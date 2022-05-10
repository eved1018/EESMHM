import pandas as pd
import plotly.express as px

def heatmap(df,pdb):
    amino_acids = ['V', 'I', 'L', 'E', 'Q', 'D', 'N', 'H', 'W', 'F', 'Y', 'R', 'K', 'S', 'T', 'M', 'A', 'G', 'P', 'C']
    cols = df.columns
    df["mutant_amino_acid"] = [i[-1] for i in df["mutant"]]
    df["wt"] = [i[:-1] for i in df["mutant"]]
    for col in cols:
        if col != "mutant":
            df2 = pd.DataFrame(index = df["wt"].unique(), columns=amino_acids)
            for index, row in df.iterrows():
                maa, wt = row["mutant"][-1], row["mutant"][:-1]
                val = row[col]
                if val != "nan":
                    df2.at[wt,maa] = float(val)

            df2 = df2.dropna(axis=1, how='all')
            fig = px.imshow(df2,text_auto = True,labels=dict(x="Mutant Amino Acid", y="Position", color=col ))
            fig.update_xaxes(side="top")
            fig.write_image(f"output/{pdb}/{col}.svg")
    return
