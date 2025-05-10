from langchain import SerpAPIWrapper
from langchain.tools import BaseTool
from langchain.llms import BaseLLM
from typing import List, Optional, Tuple
import os
from neo4j import GraphDatabase
from langchain.chains import LLMChain
from langchain.prompts import PromptTemplate

# WebSearch Tool
def web_search(keywords, search_engine="google"):
    try:
        return SerpAPIWrapper(
            serpapi_api_key=os.getenv("SERP_API_KEY"), search_engine=search_engine
        ).run(keywords)
    except:
        return "No results, try another search"

class WebSearch(BaseTool):
    name = "WebSearch"
    description = (
        "Input a specific question, returns an answer from web search. "
        "Do not mention any specific molecule names, but use more general features to formulate your questions."
    )
    serp_api_key: str = None

    def __init__(self, serp_api_key: str = None):
        super().__init__()
        self.serp_api_key = serp_api_key

    def _run(self, query: str) -> str:
        if not self.serp_api_key:
            return "No SerpAPI key found. This tool may not be used without a SerpAPI key."
        return web_search(query)

    async def _arun(self, query: str) -> str:
        raise NotImplementedError("Async not implemented")

# KnowledgeGraphSearch Tool
class KnowledgeGraphSearch(BaseTool):
    name = "KnowledgeGraphSearch"
    description = "Performs intelligent knowledge graph exploration. Input: Natural language question. Output: Relevant conclusions."
    llm: BaseLLM = None
    max_iterations: int = None
    visited_keywords: set = None

    def __init__(self, llm):
        super().__init__()
        self.llm = llm
        self.max_iterations = 5
        self.visited_keywords = set()

    def _run(self, query: str) -> str:
        """Main entry point with intelligent filtering"""
        collected_triplets = []
        current_keywords = self._generate_initial_keywords(query)
        original_query = query

        for iteration in range(1, self.max_iterations + 1):
            # Execute search and filter results
            raw_triplets = self._execute_graph_search(current_keywords)
            filtered = self._filter_relevant_triplets(raw_triplets, original_query)

            if not filtered:
                continue

            current_triplet_str = "\n".join(filtered)
            collected_triplets.append(current_triplet_str)

            # Check stopping condition
            if self._should_stop("\n".join(collected_triplets), original_query):
                break

            # Generate next keywords
            current_keywords = self._generate_next_keywords(
                original_query,
                current_triplet_str,
                iteration
            )
            self.visited_keywords.update(current_keywords)

        # Generate a structured answer from the collected triplets
        final_answer = self._summarize_triplets("\n\n".join(collected_triplets), original_query)
        return final_answer

    def _summarize_triplets(self, triplets: str, query: str) -> str:
        """
        Summarize the collected triplets into a structured answer.

        Args:
            triplets: Collected triplets as a string.
            query: The original user query.

        Returns:
            A structured answer summarizing the triplets.
        """
        prompt = PromptTemplate(
            input_variables=["triplets", "query"],
            template="""You are a knowledge graph expert. Your task is to generate a clear, logical, and structured answer based on the collected triplets to answer the user's query.

User Query: {query}

Collected Triplets:
{triplets}

Instructions:
1. Analyze the triplets and extract key information relevant to the query.
2. Organize the information into a structured format (e.g., bullet points or numbered list).
3. Ensure the answer is concise, accurate, and directly addresses the query.
4. If the triplets do not provide enough information to answer the query, state that clearly.

Output:"""
        )

        chain = LLMChain(llm=self.llm, prompt=prompt)
        response = chain.run({
            "triplets": triplets,
            "query": query
        })

        return response

    def _generate_initial_keywords(self, query: str) -> List[str]:
        """Generate context-aware initial keywords"""
        prompt = PromptTemplate(
            input_variables=["query"],
            template="""Extract 3-5 core entities from this query. 

Query: {query}

Rules:
1. Use only entity names (no verbs/adjectives).
2. Return comma-separated values.

Entities:"""
        )
        chain = LLMChain(llm=self.llm, prompt=prompt)
        response = chain.run({
            "query": query
        })
        return [kw.strip() for kw in response.split(",") if kw.strip()]

    def _generate_next_keywords(self, query: str, triplets: str, iteration: int) -> List[str]:
        """Generate refined search keywords"""
        prompt = PromptTemplate(
            input_variables=["query", "triplets", "iteration", "used_keywords"],
            template="""You are a knowledge graph expert. Current query: {query}

Collected Triplets (Iteration {iteration}):
{triplets}

Task: Generate 3-5 new search keywords to find missing links
Guidelines:
1. Focus on expanding relationships between existing entities.
2. Avoid duplicates: {used_keywords}
3. Use only comma-separated terms.

Output:"""
        )
        chain = LLMChain(llm=self.llm, prompt=prompt)
        response = chain.run({
            "query": query,
            "triplets": triplets,
            "iteration": iteration,
            "used_keywords": ", ".join(self.visited_keywords)
        })
        return [kw.strip() for kw in response.split(",") if kw.strip()]

    def _execute_graph_search(self, keywords: List[str]) -> List[str]:
        """Execute Neo4j query with parameterization"""
        if not keywords:
            return []

        triplets = []
        from neo4j import GraphDatabase
        os.environ["NEO4J_URI"] = "neo4j+s://bb60f546.databases.neo4j.io:7687"
        os.environ["NEO4J_USERNAME"] = "neo4j"
        os.environ["NEO4J_PASSWORD"] = "sWFLqNrAjD50BrArVUhQHh3CKiPSH0qJnUPU0nW1BpQ"
        driver = GraphDatabase.driver(os.getenv("NEO4J_URI"),
                                      auth=(os.getenv("NEO4J_USERNAME"), os.getenv("NEO4J_PASSWORD")))

        for keyword in keywords:
            with driver.session() as session:
                # Construct the Cypher query to search for nodes and their relationships
                cypher_query = f"""
                MATCH (n) 
                WHERE toLower(n.id) CONTAINS toLower("{keyword}")
                MATCH (n)-[r]->(m)
                RETURN DISTINCT n.id AS head, TYPE(r) AS relation, m.id AS tail
                LIMIT 100
                """
                results = session.run(cypher_query, keywords=keyword)
                # Execute the query and process the results
                for result in results:
                    triplets.append((result['head'], result['relation'], result['tail']))
        return triplets

    def _filter_relevant_triplets(self, triplets: List[Tuple[str, str, str]], query: str, top_k: int = 10) -> List[str]:
        """
        Filter the most relevant triplets using LLM.

        Args:
            triplets: List of triplets in the format [(head, relation, tail), ...].
            query: The user query.
            top_k: Number of most relevant triplets to return.

        Returns:
            List of the most relevant triplets in the format ["head - relation -> tail", ...].
        """
        if not triplets:
            return []

        # Convert triplets to string format
        triplet_strs = [f"{h} - {r} -> {t}" for h, r, t in triplets]

        # Construct prompt template
        prompt = PromptTemplate(
            input_variables=["triplets", "query", "top_k"],
            template="""Evaluate the relevance of the following triplets to the query and select the top {top_k} most relevant ones.

Query: {query}

Triplets:
{triplets}

Instructions:
1. Consider direct relevance to query entities.
2. Consider logical connection to the query intent.
3. Consider information value for answering the query.

Output format:
- Return only the selected triplets, one per line.
- Do not include any additional text or explanations.
"""
        )

        # Construct input
        chain = LLMChain(llm=self.llm, prompt=prompt)
        response = chain.run({
            "triplets": "\n".join(triplet_strs),
            "query": query,
            "top_k": top_k
        })

        # Parse response
        selected_triplets = [line.strip() for line in response.strip().split("\n") if line.strip()]
        return selected_triplets[:top_k]  # Ensure no more than top_k triplets

    def _should_stop(self, all_triplets: str, query: str) -> bool:
        """LLM-based stopping condition"""
        if not all_triplets:
            return False

        prompt = PromptTemplate(
            input_variables=["triplets", "query"],
            template="""Determine if sufficient to answer:

Query: {query}
Collected Triplets:
{triplets}

Analysis Factors:
1. Key entity coverage
2. Relationship completeness
3. Missing information

Decision (one word): Sufficient/Insufficient"""
        )

        chain = LLMChain(llm=self.llm, prompt=prompt)
        response = chain.run({
            "triplets": all_triplets,
            "query": query
        }).strip()
        return response == "Sufficient"

    async def _arun(self, query: str) -> str:
        raise NotImplementedError("Async operation not supported")