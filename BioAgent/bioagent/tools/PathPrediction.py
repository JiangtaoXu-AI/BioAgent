"""Self-hosted reaction tools. Retrosynthesis, reaction forward prediction."""

from typing import Union

from rdkit import Chem

from typing import List, Dict, Any

from typing import Optional
from langchain.llms import BaseLLM

from langchain import LLMChain, PromptTemplate

from .prompts import prediction_info_prompt,final_info_prompt

from langchain.tools import BaseTool
import requests

from pydantic import Field
from langchain.memory import ConversationBufferMemory

__all__ = ["SMILESToPredictedSynthesisInfo","AddReactantsToBionavi","SMILESToBiosyntheticPathway"]
class SMILESToPredictedSynthesisInfo(BaseTool):
    """
    SMILES 转合成相关信息工具（内置 LLM 提取关键信息）
    """
    name: str = "SMILESToPredictedSynthesisInfo"
    description: str = "Input target molecule SMILES, and return synthesis related information through the Bionavi library, including Enzymes, reactants, reactions, and pathways."

    llm: BaseLLM = None
    llm_chain1: LLMChain = None
    llm_chain2: LLMChain = None
    # 显式声明 memory 字段
    memory: Optional[ConversationBufferMemory] = Field(default=None)

    def __init__(self, llm: BaseLLM,memory):
        super().__init__()
        self.llm = llm
        self.memory = memory  # 共享内存
        prompt1 = PromptTemplate(
            template=prediction_info_prompt,
            input_variables=["data"]
        )
        self.llm_chain1 = LLMChain(prompt=prompt1, llm=self.llm)
        prompt2 = PromptTemplate(
            template=final_info_prompt,
            input_variables=["data", "stored_data"]
        )
        self.llm_chain2 = LLMChain(prompt=prompt2, llm=self.llm)

    def predict_synthesis_info(self, smiles: str) -> str:
        # FastAPI 服务的 URL
        url = "http://101.33.233.196/api/v1/plan"

        # 请求体
        payload = {"smiles": smiles}

        try:
            # 发送 POST 请求
            response = requests.post(url, json=payload)
            response.raise_for_status()  # 检查请求是否成功

            # 解析 API 响应
            result = response.json()

            # 格式化输出
            output = []
            for route in result:
                output.append(f"route id: {route['route_id']}")
                output.append(f"route score: {route['route_score']}")
                output.append(f"route: {route['route']}")
                output.append("")  # 添加空行分隔

            return "\n".join(output)

        except requests.exceptions.RequestException as e:
            return f"Error: {str(e)}"

    def _run(self, smiles: str) -> str:
        """
        输入 SMILES，返回提取后的合成相关信息
        :param smiles: 化合物的 SMILES 表示
        :return: 提取后的合成相关信息或提示信息（字符串）
        """
        try:
            # 调用预测函数获取原始输出
            print((f"\nPredicting the synthesis pathway of {smiles}. This will take approximately 1 minute, please be patient........"))
            raw_output = self.predict_synthesis_info(smiles)

            # 尝试从内存中获取最近一次存储的数据
            stored_data = ""
            if self.memory and self.memory.buffer:
                last_entry = self.memory.buffer_as_messages[-1].content
                if isinstance(last_entry, dict):  # 确保是字典
                    stored_data = last_entry.get("output", "")
                elif isinstance(last_entry, str):  # 如果是字符串，直接使用
                    stored_data = last_entry

            # 如果没有返回任何数据，则返回统一提示信息
            if not raw_output.strip():
                return "No synthesis information found for the given molecule."+"Final synthesized information\n"+stored_data
            # 使用 LLM 解析并提取关键信息
            #print(f"stored_data :{stored_data}")
            prediction_output = self.llm_chain1.run(data=raw_output)
            #print(f"prediction_output:{prediction_output}")
            final_output = self.llm_chain2.run(data=prediction_output,stored_data=stored_data)
            return final_output
        except Exception as e:
            return f"An error occurred while processing the request: {str(e)}"


class AddReactantsToBionavi(BaseTool):
    """
    将反应物添加到 Bionavi 规则库的工具
    """
    name: str = "AddReactantsToBionavi"
    description: str = (
        "Before using the SMILESToPredictedSynthesisInfo tool to predict synthesis routes for target molecules, "
        "this tool allows you to enhance the accuracy of predictions by adding potential reactants to the Bionavi rule library. "
        "The input is in the following form: Organisms=[...],Enzymes=[...],Reactants=[...], Reactions=[...], Pathways=[...](e.g., Organisms=[Archaea,.........])."
    )
    # 显式声明 memory 字段
    memory: Optional[ConversationBufferMemory] = Field(default=None)

    def __init__(self,memory):
        super().__init__()
        self.memory=memory
    def _run(self, input_data) -> str:
        """
        输入是一个逗号分割的字符串，处理并添加到规则库
        :param input_data: 逗号分割的分子字符串 (e.g., 'O=C=O, [CH3]')
        :return: 添加结果的描述信息
        """
        if not input_data or not isinstance(input_data, str):
            return "Invalid input. Please provide a comma-separated string of reactants."

        # 分割字符串为列表，并清理每个元素
        reactants_list = [r.strip() for r in input_data.split(",") if r.strip()]

        # 假设这里执行了将反应物添加到规则库的操作
        num_added = len(reactants_list)

        # 存储数据到内存
        if self.memory:
            self.memory.save_context({"input": input_data}, {"output": input_data})

        # 返回具体的结果信息
        return f"Successfully added {num_added} potential reactants to the Bionavi rule library."

    async def _arun(self, input_data) -> str:
        """
        异步版本的运行方法（尚未实现异步逻辑）
        """
        raise NotImplementedError("Async version not implemented.")
class SMILESToBiosyntheticPathway(BaseTool):
    """
    SMILES 转 Biosynthetic Pathway 工具
    """
    name: str = "SMILESToBiosyntheticPathway"
    description: str = "Input a molecule's SMILES, returns the biosynthetic pathway in the current database."
    url: str = None

    def __init__(self):
        super().__init__()
        self.url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{}/cids/JSON"

    def validate_smiles(self, smiles: str) -> bool:
        """
        验证 SMILES 格式是否正确
        :param smiles: 化合物的 SMILES 表示
        :return: 是否合法
        """
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None

    def get_compound_id_from_smiles(self, smiles: str) -> Union[str, None]:
        """
        通过 SMILES 查询 PubChem 获取化合物 ID (CID)
        :param smiles: 化合物的 SMILES 表示
        :return: 化合物 CID 或 None
        """
        url = self.url.format(smiles)
        try:
            response = requests.get(url, timeout=10)
            if response.status_code == 200:
                data = response.json()
                cids = data.get("IdentifierList", {}).get("CID", [])
                if cids:
                    return str(cids[0])
        except Exception:
            pass
        return None

    def get_kegg_compound_id_from_pubchem(self, cid: str) -> Union[str, None]:
        """
        通过 PubChem CID 查询对应的 KEGG 化合物 ID
        :param cid: PubChem 化合物 ID
        :return: KEGG 化合物 ID 或 None
        """
        url = f"http://rest.kegg.jp/conv/compound/pubchem:{cid}"
        try:
            response = requests.get(url, timeout=10)
            if response.status_code == 200:
                lines = response.text.split("\n")
                for line in lines:
                    if line.startswith("cpd:"):
                        return line.split(":")[1]
        except Exception:
            pass
        return None

    def find_pathway_with_compound(self, compound_id: str) -> Union[str, None]:
        """
        查找包含指定化合物的代谢通路
        :param compound_id: KEGG 化合物 ID
        :return: 通路 ID 或 None
        """
        url = f"http://rest.kegg.jp/link/pathway/{compound_id}"
        try:
            response = requests.get(url, timeout=10)
            if response.status_code == 200:
                lines = response.text.split("\n")
                for line in lines:
                    if line.startswith("path:"):
                        return line.split(":")[1]
        except Exception:
            pass
        return None

    def generate_reaction_chain(self, pathway_id: str) -> Union[List[Dict[str, Any]], None]:
        """
        从指定通路中生成反应链条
        :param pathway_id: KEGG 通路 ID
        :return: 反应链条列表 或 None
        """
        url = f"http://rest.kegg.jp/get/{pathway_id}"
        try:
            response = requests.get(url, timeout=10)
            if response.status_code != 200:
                return None

            reactions = []
            lines = response.text.split("\n")
            in_reaction_section = False
            for line in lines:
                if line.startswith("REACTION"):
                    in_reaction_section = True
                    continue
                if not line.strip() or not in_reaction_section:
                    continue
                if line.startswith(" "):  # 反应信息通常缩进
                    reaction_ids = line.strip().split()
                    reactions.extend(reaction_ids)
                else:
                    break  # 结束反应部分

            chain = []
            for reaction_id in reactions:
                reaction_details = self.get_reaction_details(reaction_id)
                if reaction_details:
                    chain.append(reaction_details)

            return chain if chain else None
        except Exception:
            return None

    def get_reaction_details(self, reaction_id: str) -> Union[Dict[str, Any], None]:
        """
        获取单个反应的详细信息
        :param reaction_id: 反应 ID
        :return: 反应详情字典 或 None
        """
        url = f"http://rest.kegg.jp/get/{reaction_id}"
        try:
            response = requests.get(url, timeout=10)
            if response.status_code != 200:
                return None

            details = {"id": reaction_id, "equation": ""}
            lines = response.text.split("\n")
            in_equation_section = False
            for line in lines:
                if line.startswith("EQUATION"):
                    in_equation_section = True
                    details["equation"] = line[8:].strip()
                    continue
                if in_equation_section and line.startswith(" "):
                    details["equation"] += line.strip()

            return details if details["equation"] else None
        except Exception:
            return None

    def _run(self, smiles: str) -> str:
        """
        输入 SMILES，返回可能的生物合成路径或提示信息
        :param smiles: 化合物的 SMILES 表示
        :return: 反应路径或提示信息（字符串）
        """
        try:
            if not self.validate_smiles(smiles):
                return "Biosynthetic pathway not found in the database for this compound."

            # Step 1: 获取 PubChem CID
            cid = self.get_compound_id_from_smiles(smiles)
            if not cid:
                return "Biosynthetic pathway not found in the database for this compound."

            # Step 2: 获取 KEGG 化合物 ID
            kegg_compound_id = self.get_kegg_compound_id_from_pubchem(cid)
            if not kegg_compound_id:
                return "Biosynthetic pathway not found in the database for this compound."

            # Step 3: 查找包含该化合物的代谢通路
            pathway_id = self.find_pathway_with_compound(kegg_compound_id)
            if not pathway_id:
                return "Biosynthetic pathway not found in the database for this compound."

            # Step 4: 生成反应链条
            reaction_chain = self.generate_reaction_chain(pathway_id)
            if not reaction_chain:
                return "Biosynthetic pathway not found in the database for this compound."

            # 如果找到路径，将反应链条格式化为字符串输出
            pathway_str = "\n".join(
                [f"Reaction {i}: {reaction['id']}\nEquation: {reaction['equation']}"
                 for i, reaction in enumerate(reaction_chain, 1)]
            )
            return pathway_str
        except Exception:
            return "Biosynthetic pathway not found in the database for this compound."



