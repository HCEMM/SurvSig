<p align="center">
  <img src="https://github.com/user-attachments/assets/1ac818d1-0cbe-49d3-ad4f-764d3add6c5b" width="500">
</p>

# **SurvSig**  
## *A web-based tool for complex gene signature analysis in neuroendocrine lung tumors and TCGA*

📞 **Access the application here:** [**www.survsig.hcemm.eu**](https://survsig.hcemm.eu/)  
💡 **No installation required** – just open the website and start analyzing!  
🔬 **About our group:** [**www.hcemm.eu**](https://www.hcemm.eu/teams/genomic-instability-and-cancer/cancer-genomics-and-epigenetics-core-group/)

---

## 🚀 **About SurvSig**  
SurvSig is an **interactive web application** designed for **clinicians and researchers** working with **neuroendocrine lung tumors** (*SCLC, LCNEC, carcinoid tumors*).  
It enables users to analyze **complex gene signatures** and explore their **relationship with patient outcomes** through intuitive visualizations.  

### 🛡️ **Why use SurvSig?**  
✅ **Designed for clinical research** – No coding or bioinformatics expertise required.  
✅ **Supports multiple real-world datasets** – Including **SCLC, LCNEC, carcinoid tumors, and TCGA**.  
✅ **User-friendly, web-based interface** – Simply open the website and start analyzing.  
✅ **Gene signature-based analysis** – Upload your own gene list or use predefined signatures.  
✅ **Advanced visualization tools** – Explore survival plots, heatmaps, and UMAP clustering.  

---

## 📊 **Key Features**  
✔ **Multi-cohort analysis** – Compare multiple patient datasets in one platform.  
✔ **Support for major datasets:**  
   - **🦰 Small Cell Lung Cancer (SCLC):** *George-SCLC, Liu-SCLC, Lissa-SCLC, Jiang-SCLC*  
   - **🧬 Large-Cell Neuroendocrine Carcinoma (LCNEC):** *George-LCNEC*  
   - **🩺 Carcinoid Tumors:** *Alcala-Carcinoid, Fernandez-Carcinoid*  
   - **🌐 Mixed Cohort Data:** *Rousseaux Lung Tumors, TCGA*  
✔ **Machine learning-powered insights** – Identify molecular subtypes with **advanced clustering methods**.  
✔ **Custom gene list support** – Use **predefined or user-defined** gene signatures for analysis.  

<p align="center">
  <img src="https://github.com/user-attachments/assets/f684b2aa-9bcf-4ae7-a13c-24766a18db9f" width="400">
</p>

---

## 🛠️ **How to Use SurvSig (Online)**  
1️⃣ **Visit the website:** [www.survsig.hcemm.eu](https://survsig.hcemm.eu/)  
2️⃣ **Select a dataset and analysis type.**  
3️⃣ **Upload a gene list** or choose from predefined signatures.  
4️⃣ **Explore results**: View **survival analysis, enrichment scores, heatmaps, and clustering outputs.**  

---

## 🖥 **How to Install Locally (Optional)**  
If you want to run SurvSig locally and upload custom datasets, follow these steps:

### **Prerequisites**
- Python **>= 3.8**  
- Git installed (`git --version` to check)  
- Pipenv (recommended) or pip  

### **Installation Steps**
1️⃣ **Clone the repository**  
```bash
git clone https://github.com/HCEMM/SurvSig.git
cd SurvSig
```

2️⃣ **Set up the environment**  
Using `pipenv` (recommended):  
```bash
pip install pipenv
pipenv install
pipenv shell
```
Or using `pip`:  
```bash
pip install -r requirements.txt
```

3️⃣ **Run the application**  
```bash
streamlit run main.py
```

4️⃣ **Open the application**  
- After running the command, the app will start at:  
  **`http://localhost:8501/`**  

### **How to Add Custom Data?**
- Place your **gene expression matrix** and **clinical data** into the `/source_data` folder.  
- Integrate your dataset into the code

---

## 📚 **Citation**  
If you use **SurvSig** in your research, please cite:  
📚 *Nemes et al., 2025 (manuscript in preparation)*  

---

## 📝 **License**  
SurvSig is developed at **HCEMM** and is available under the **GLP v3**.  
