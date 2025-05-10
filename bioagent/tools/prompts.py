safety_summary_prompt = (
    "Your task is to parse through the data provided and provide a summary of important health, laboratory, and environemntal safety information."
    'Focus on answering the following points, and follow the format "Name: description".'
    "Operator safety: Does this substance represent any danger to the person handling it? What are the risks? What precautions should be taken when handling this substance?"
    "GHS information: What are the GHS signal (hazard level: dangerous, warning, etc.) and GHS classification? What do these GHS classifications mean when dealing with this substance?"
    "Environmental risks: What are the environmental impacts of handling this substance."
    "Societal impact: What are the societal concerns of this substance? For instance, is it a known chemical weapon, is it illegal, or is it a controlled substance for any reason?"
    "For each point, use maximum two sentences. Use only the information provided in the paragraph below."
    "If there is not enough information in a category, you may fill in with your knowledge, but explicitly state so."
    "Here is the information:{data}"
)

summary_each_data = (
    "Please summarize the following, highlighting important information for health, laboratory and environemntal safety."
    "Do not exceed {approx_length} characters. The data is: {data}"
)
prediction_info_prompt= (
    "Input format:"
    "Predictive Information: {data}\n"  # 使用 {data} 占位符
    "Your task is to parse through the Predictive Information provided and extract synthesis-related information."
    "The Predictive Information contains multiple routes with IDs, scores, and detailed steps"
    "Summarize all routes of the Predictive Information into 4 aspects: enzymes, reactants, reactions and pathways according to the following format:\n"
    "For each aspect(enzymes,reactants,reactions,pathways) first provide a summary sentence, followed by the components listed as comma-separated values.\n\n"
    "Example format for Reactions: \n"
    "Reactions: Mainly phosphorylation and dehydrogenation reactions; [C:1]-[OH;D1;+0:2] >> O-P(-O)(=O)-[O;H0;D2;+0:2]-[C:1], [C:1]-[CH2;D2;+0:2]-[CH;D3;+0:3](-[C:4])-[C:5] >> [C:1]-[CH;D2;+0:2]=[C;H0;D3;+0:3](-[C:4])-[C:5]...\n"
    "Example format for pathways: \n"
    "route0: Key Overview of Pathway; smile1>smile2>smile3........\n route1:................."
)


final_info_prompt = (
    "**Input:"
    "Predictive Information: {data}\n"
    "Literature information: {stored_data}\n"
    "**Your task is:"
    "You are a biochemist tasked with integrating both Literature Information and Predictive Information to obtain more accurate and complete synthesis information.\n\n"
    "The Literature Information includes Organisms, Enzymes, Reactants, Reactions, and Pathways, which are mainly in textual form and are considered reliable.\n"
    "The Predictive Information includes Enzymes, Reactants, Reactions, and Pathways. For each aspect, contains rich symbolic information such as the smile structure\n\n"
    "You need to combine both sources of information, Considering the relationships between textual and symbolic information and biochemical knowledge, to generate more accurate and complete synthesis data.\n\n"
    "**Output format:\n"
    "5 aspects: Organisms, Enzymes, Reactants, Reactions, Pathways.\n"
    "For each aspect, first provide a summary sentence, followed by the components in a comma-separated list (including both textual information and symbolic notation).\n"
    "**Sample output for enzyme:\n"
    "Enzymes:\n"
    "Summary: Mainly related to the biosynthesis of macrolide compounds, lipids, and secondary metabolites.\n"
    "Components: branched alkyltransferase, GGGPS, PlsC, MACROLIDE-2-KINASE, RXN18C3-240, MetaCyc_MACROLIDE-2-KINASE-RXN, MetaCyc_RXN18C3-240\n"
    "Pathways:"
    "Summary: The main idea behind these paths is......"
    "route0: Key Overview of route; smile1>smile2>smile3........\n "
    "route1:................."
)