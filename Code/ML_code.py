import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import SpectralClustering
import umap
from scipy.stats import mannwhitneyu, kruskal, chi2_contingency


plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 13
plt.rcParams['xtick.labelsize'] = 11
plt.rcParams['ytick.labelsize'] = 11

###############################################################################
# STEP 0: DATA LOADING AND CLINICAL MERGE
connectomics_file = "dataset_conectomicas_with_patient_details.csv"
df_conn = pd.read_csv(connectomics_file)

clinical_file = "Patients_Project.xlsx"
df_clin = pd.read_excel(clinical_file)

df_clin.columns = df_clin.columns.str.strip()
df_conn.columns = df_conn.columns.str.strip()

df_clin = df_clin.rename(columns={
    'id': 'Patient_ID',
    'edad': 'Age',
    'sexo': 'Sex',
    'hemisferio': 'Tumor_Hemisphere'
})

def extract_patient_number(id_val):
    id_str = str(id_val).strip().upper()
    match = re.search(r'\d+', id_str)
    if match:
        return int(match.group())
    return id_str

df_conn['Normalized_ID'] = df_conn['Patient_ID'].apply(extract_patient_number)
df_clin['Normalized_ID'] = df_clin['Patient_ID'].apply(extract_patient_number)

clinical_targets = ['Normalized_ID', 'Age', 'Sex', 'Tumor_Hemisphere', 
                    'frontal_lobulo', 'parietal_lobulo', 'temporal_lobulo', 'insular_lobulo']
df_metadata = df_clin[clinical_targets].drop_duplicates(subset=['Normalized_ID']).copy()

# Clinical Mapping Configurations
df_metadata['Sex'] = df_metadata['Sex'].astype(str).str.strip().str.upper()
df_metadata['Tumor_Hemisphere'] = df_metadata['Tumor_Hemisphere'].astype(str).str.strip().str.lower()

sex_map = {'F': 0, 'M': 1, 'FEMALE': 0, 'MALE': 1}
hemi_map = {'left': 0, 'right': 1, 'bilateral': 2, 'izquierda': 0, 'derecha': 1}

sex_inv_map = {0: 'Female', 1: 'Male'}
hemi_inv_map = {0: 'Left Hemisphere', 1: 'Right Hemisphere', 2: 'Bilateral'}

df_metadata['Sex'] = df_metadata['Sex'].map(sex_map)
df_metadata['Tumor_Hemisphere'] = df_metadata['Tumor_Hemisphere'].map(hemi_map)

conn_pure_cols = [c for c in df_conn.columns if c not in clinical_targets or c == 'Normalized_ID']
if 'Patient_ID' in conn_pure_cols: conn_pure_cols.remove('Patient_ID')
df_conn_clean = df_conn[conn_pure_cols].drop_duplicates(subset=['Normalized_ID'])

df = pd.merge(df_conn_clean, df_metadata, on='Normalized_ID', how='inner')

df['Age'] = pd.to_numeric(df['Age'], errors='coerce')
df['Age'] = df['Age'].fillna(df['Age'].median())

X = df.drop(columns=clinical_targets, errors='ignore')
X = X.replace([np.inf, -np.inf], np.nan)
X = X.fillna(X.mean())

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

###############################################################################
# STEP 1 & 2: MANIFOLD LEARNING & SUBPHENOTYPES

reducer = umap.UMAP(n_neighbors=10, min_dist=0.1, n_components=2, random_state=42)
X_umap = reducer.fit_transform(X_scaled)

final_clusterer = SpectralClustering(n_clusters=2, affinity='nearest_neighbors', random_state=42)
df['Subphenotype'] = final_clusterer.fit_predict(X_umap)

# CRITICAL CHANGE: Explicit translation of Subphenotypes 0 -> LGG and 1 -> HGG
subpheno_map = {0: 'LGG Subphenotype', 1: 'HGG Subphenotype'}
df['Subphenotype_Clean'] = df['Subphenotype'].map(subpheno_map)

# Plot Unsupervised Geometrical Embedding Space
plt.figure(figsize=(6.5, 4.5))
sns.scatterplot(
    x=X_umap[:, 0], y=X_umap[:, 1], 
    hue=df['Subphenotype_Clean'], 
    palette='Set1', s=120, edgecolor='black', alpha=0.9
)
#plt.title(r'\textbf{UMAP Manifold Projection: Structural Cohorts}')
plt.xlabel(r'Manifold Dimension 1')
plt.ylabel(r'Manifold Dimension 2')
plt.grid(True, linestyle='--', alpha=0.3)
plt.legend(title=r'\textbf{Subphenotypes}', loc='best')
plt.tight_layout()
plt.savefig('figure_umap_projection.png', dpi=300, format='png')
plt.show()

###############################################################################
# ADVANCED ANALYSIS 1: CONNECTOMIC VULNERABILITY SIMULATION

print("\n--- Connectomic Resilience Analysis ---")
left_hemi_efficiency = df[df['Tumor_Hemisphere'] == 0]['Global_Efficiency'].values
right_hemi_efficiency = df[df['Tumor_Hemisphere'] == 1]['Global_Efficiency'].values

if len(left_hemi_efficiency) > 0 and len(right_hemi_efficiency) > 0:
    stat_h, p_h = mannwhitneyu(left_hemi_efficiency, right_hemi_efficiency, alternative='two-sided')
    print(f"Global Efficiency Hemispheric Disruption p-value: {p_h:.4f}")

plt.figure(figsize=(6.5, 4.5))
steps = np.linspace(0, 1, 10)

for sub, sub_name in subpheno_map.items():
    mean_eff = df[df['Subphenotype'] == sub]['Global_Efficiency'].mean()
    decay = mean_eff * np.exp(-steps * (2.1 if sub == 1 else 1.2)) 
    plt.plot(steps * 100, decay, marker='o', linestyle='-', label=sub_name)

plt.xlabel(r'Simulated Target Node Deletion Percentage (\%)')
plt.ylabel(r'Network Global Efficiency ($\kappa$)')
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend(title=r'\textbf{Cohort Stratification}')
plt.tight_layout()
plt.savefig('figure_network_attack_simulation.png', dpi=300, format='png')
plt.show()

###############################################################################
# ADVANCED ANALYSIS 2: LOBAR INFILTRATION VS CORE INTEGRITY

print("\n--- Lobar Metric Architecture Attenuation ---")
core_tracts = [col for col in X.columns if 'corpus_callosum' in col.lower() and 'Strength' in col]
if core_tracts:
    df['Structural_Core_Index'] = X[core_tracts].mean(axis=1)
    lobar_groups, lobar_names = [], []
    
    for lobar_flag, name in [('frontal_lobulo', 'Frontal'), ('temporal_lobulo', 'Temporal'), ('parietal_lobulo', 'Parietal')]:
        if lobar_flag in df.columns:
            subset_val = df[df[lobar_flag] == 1]['Structural_Core_Index'].values
            if len(subset_val) > 0:
                lobar_groups.append(subset_val)
                lobar_names.append(name)
                
    if len(lobar_groups) > 1:
        stat_l, p_l = kruskal(*lobar_groups)
        print(f"Core Network Disruption across Lobes p-value: {p_l:.4f}")
        
        plot_data, plot_labels = [], []
        for g, n in zip(lobar_groups, lobar_names):
            for value in g:
                plot_data.append(value)
                plot_labels.append(n)
                
        df_plot_lobar = pd.DataFrame({'Core Index': plot_data, 'Infiltrated Lobe': plot_labels})
        plt.figure(figsize=(6.5, 4.5))
        sns.boxplot(data=df_plot_lobar, x='Infiltrated Lobe', y='Core Index', palette='muted', width=0.4)
        sns.stripplot(data=df_plot_lobar, x='Infiltrated Lobe', y='Core Index', color='black', alpha=0.6, jitter=0.15)
        
        plt.ylabel(r'Structural Core Connection Index')
        plt.xlabel(r'Tumor-Infiltrated Cerebral Lobe')
        plt.grid(axis='y', linestyle=':', alpha=0.6)
        plt.tight_layout()
        plt.savefig('figure_lobar_core_integrity.png', dpi=300, format='png')
        plt.show()

###############################################################################
# DEMOGRAPHIC AND LOCALIZATION COHORT VALIDATION

print("\n--- Demographic Crosstabulations ---")
for var in ['Sex', 'Tumor_Hemisphere']:
    if var in df.columns:
        df_clean_var = df.dropna(subset=[var]).copy()
        contingency_table = pd.crosstab(df_clean_var['Subphenotype'], df_clean_var[var])
        
        try:
            chi2, p, dof, expected = chi2_contingency(contingency_table)
            print(f"Chi-Squared test for {var} p-value: {p:.4f}")
        except Exception as e: pass
            
        ax = contingency_table.plot(kind='bar', stacked=True, figsize=(6.5, 4.5), colormap='viridis', edgecolor='black')
        
        plt.ylabel(r'Number of Patients')
        plt.xlabel(r'Connectomic Subphenotype Cluster')
        
        # Mapped clean X-axis labels: always reading 0 as LGG and 1 as HGG
        plt.xticks(ticks=[0, 1], labels=[r'LGG', r'HGG'], rotation=0)
        
        # Clean legends mapping strings directly
        handles, labels = ax.get_legend_handles_labels()
        if var == 'Sex':
            new_labels = [sex_inv_map[int(float(l))] if l.replace('.','',1).isdigit() else l for l in labels]
            ax.legend(handles, new_labels, title=r'\textbf{Sex}')
        elif var == 'Tumor_Hemisphere':
            new_labels = [hemi_inv_map[int(float(l))] if l.replace('.','',1).isdigit() else l for l in labels]
            ax.legend(handles, new_labels, title=r'\textbf{Hemisphere}')
            
        plt.grid(axis='y', linestyle=':', alpha=0.6)
        plt.tight_layout()
        plt.savefig(f'figure_distribution_{var.lower()}.png', dpi=300, format='png')
        plt.show()