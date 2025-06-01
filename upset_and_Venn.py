import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from upsetplot import UpSet, from_indicators
from venn import venn

# File and output paths
matrix_file = "/home/shared2/misbah/all_strains/rna_seq/venn/Mn4/venn_matrix.csv"
upset_out = "/home/shared2/misbah/all_strains/rna_seq/venn/Mn4/upset_plot_Mn4.png"
venn_out = "/home/shared2/misbah/all_strains/rna_seq/venn/Mn4/venn4_colored_Mn4.png"

# 1. Load Data
df = pd.read_csv(matrix_file)
indicators = df.drop(columns=['Accession'])
indicators_bool = indicators.astype(bool)

# 2. UpSet Plot (No element_color)
plt.figure(figsize=(10, 6))
upset = UpSet(
    from_indicators(indicators_bool.columns, indicators_bool),
    show_counts=True,
    subset_size='count'
)
upset.plot()
plt.suptitle('UpSet Plot: Overlapping Significant Genes', fontsize=16)
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig(upset_out)
plt.close()
print(f"Saved UpSet plot as {upset_out}")

# 3. 4-way Venn Diagram (Colorful)
set_dict = {col: set(df.loc[df[col] == 1, 'Accession']) for col in indicators.columns}
plt.figure(figsize=(10,10))
venn(set_dict, cmap=plt.cm.Set3)
plt.title("Venn Diagram of Differentially Expressed Genes", fontsize=16)
plt.tight_layout()
plt.savefig(venn_out)
plt.close()
print(f"Saved 4-way Venn diagram as {venn_out}")
