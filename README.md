# Pipeline-PGPR_AgroBioTech

#### A sequência de scripts a serem compilados neste pipeline é: 

1. download_e_organizacao_genomas.py
2. blastp.py
3. analisar_blastp.py
4. comparacao_genomas.py
5. construcao_banco_AgroBioTech.py

Para compilar o script 1, recomenda-se utilizar um cluster. Nesta pesquisa, foi utilizado os equipamentos do Laboratório Multiusuários CCCT-CP, que utiliza o software SLURM para o gerenciamento de filas de trabalho. Utilizamos o script shell "script_python_quali.sh" para compilar o script 1.

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------


#### Bibliotecas Python e softwares necessários:

1. Biblioteca Biopython - https://biopython.org/wiki/Download
2. BLAST+ - https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------


### Fluxograma do Pipeline 

<div align="center">
<img src="https://user-images.githubusercontent.com/102994978/173956756-549af94a-af76-4c41-9f99-bbabe8d9c354.png" width="600px" /> 
</div>
