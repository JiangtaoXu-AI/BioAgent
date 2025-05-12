# 一、BioAgent - 快速入门
## 1. 克隆仓库
···python
git clone https://github.com/JiangtaoXu-AI/BioAgent.git
cd BioAgent

## 2. 创建虚拟环境
conda create -n bioagent3.10 python=3.10 -y
conda activate bioagent3.10
pip install -r requirements.txt


## 3. 配置 API 密钥
os.environ["OPENAI_API_KEY"] = "你的API密钥"
os.environ["OPENAI_API_BASE"] = "https://api.deepseek.com/v1"
model = "deepseek-chat"

## 4.运行
python main.py

# 二、工具添加
## 第一步：工具实现

- 工具的实现应放在 `bioagent/tools` 目录下。
- 每一类工具可以单独创建一个 `.py` 文件，例如 `my_tool.py`。
- 模仿别的工具格式，在该文件内编写工具的核心逻辑。

## 第二步：工具导入

在 `bioagent/agents/tools` 目录下，将工具导入并注册到工具集合中。
all_tools += [
    KnowledgeGraphSearch(llm=llm),
    SMILESToBiosyntheticPathway(),
    AddReactantsToBionavi(memory),
    SMILESToPredictedSynthesisInfo(llm=llm,memory=memory),
    GenomeCollectorTool(),
    GenomeQueryTool()
]

# 三、工具测试

## 方法一：在实现文件内测试

- 可以直接在实现文件中添加 `__main__` 代码块进行简单测试。
- 适合本地调试和单元测试使用。
if __name__ == "__main__":
    result = example_tool("apple")
    print(result)  # 输出：Processed apple
## 方法二：在 Agent 系统中测试
- 启动 agent对话框（参考1快速开始python main.py）
- 在对话框输入触发工具的用户需求
- 测试结果




