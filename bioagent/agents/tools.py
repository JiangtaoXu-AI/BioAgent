from langchain import agents
from langchain.base_language import BaseLanguageModel
from bioagent.tools import *
from langchain.memory import ConversationBufferMemory

def make_tools(llm: BaseLanguageModel, api_keys: dict = {}, local_rxn: bool=False, verbose=True):


    # 创建共享内存
    memory = ConversationBufferMemory(memory_key="chat_history", return_messages=False)
    chemspace_api_key = api_keys.get("CHEMSPACE_API_KEY") or os.getenv(
        "CHEMSPACE_API_KEY"
    )
    all_tools = agents.load_tools([])
    all_tools += [
        KnowledgeGraphSearch(llm=llm),
        SMILESToBiosyntheticPathway(),
        AddReactantsToBionavi(memory),
        SMILESToPredictedSynthesisInfo(llm=llm,memory=memory),
        GenomeCollectorTool(),
        GenomeQueryTool(),
        SelectRepresentativeSeq(),
    ]

    return all_tools
