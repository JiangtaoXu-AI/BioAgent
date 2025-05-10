from bioagent.agents import BioAgent
import os
os.environ["OPENAI_API_KEY"] = ""
os.environ["OPENAI_API_BASE"] = "https://api.deepseek.com/v1"
model="deepseek-chat"


# 初始化agent
chem_model = BioAgent(model=model, tools_model=model, temp=0.1)

KG_information=""
# 对话历史存储
history = []

while True:
    try:
        # 获取用户输入
        user_input = input("\nUser: ")

        if user_input.lower() in ["exit", "quit"]:
            break

        # 运行agent
        response = chem_model.run(user_input, history=history)

        # 显示响应并保存历史
        print(f"\nAgent: {response}")
        history.append((user_input, response))

    except KeyboardInterrupt:
        break


