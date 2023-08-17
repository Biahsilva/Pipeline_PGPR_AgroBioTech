from Bio import Entrez, SeqIO
import sys
import os

Entrez.email = 'beatrizsilva@alunos.utfpr.edu.br'
Entrez.api_key = '3d817b61bb5164b63368b92a9f158f8d4308'

# Função para ler o arquivo Input
def ler_arquivo(nome_arquivo):
    try:
        with open(nome_arquivo, 'r') as arquivo:
            linhas = arquivo.readlines()
            valores = [linha.strip() for linha in linhas]
            return valores
    except FileNotFoundError:
        print("Arquivo não encontrado.")
        return []

# Função para criar todas as pastas necessárias para o funcionamento do Pipeline
def criar_pastas():

    os.makedirs('../Banco_de_dados')
    os.makedirs('../Banco_de_dados/banco_de_dados_genomas')
    os.makedirs('../Banco_de_dados/banco_de_dados_adicionais')
    os.makedirs('../Banco_de_dados/banco_de_dados_genes/banco_de_dados_geral')
    os.makedirs('../Banco_de_dados/banco_de_dados_genes/banco_de_dados_pgp')
    os.makedirs('../Banco_de_dados/banco_de_dados_blast/banco_de_dados_sbjcts')
    os.makedirs('../Banco_de_dados/banco_de_dados_blast/banco_de_dados_queries')
    os.makedirs('../Banco_de_dados/banco_de_dados_AgroBioTech')
    os.makedirs('../Banco_de_dados/banco_de_dados_AgroBioTech/genes_pgp_referencia')
    os.makedirs('../Banco_de_dados/banco_de_dados_AgroBioTech/genes_pgp_adicionais')
    os.makedirs('../Banco_de_dados/banco_de_dados_AgroBioTech/banco_de_genes_AgroBioTech')
    os.makedirs('../Relatorios_blast')
    os.makedirs('../Relatorios_blast/relatorios_blast_complementar')
    os.makedirs('../Relatorios_blast/relatorios_blast_comparacao')
    os.makedirs('../Relatorios_blast/relatorios_blast_comparacao/90')
    os.makedirs('../Relatorios_blast/relatorios_blast_comparacao/80')
    os.makedirs('../Relatorios_blast/relatorios_blast_comparacao/70')
    os.makedirs('../Relatorios_blast/relatorios_blast_comparacao/60')
    os.makedirs('../Resultados_comparacao')
    os.makedirs('../Resultados_similaridade')
    os.makedirs('../Resultados_similaridade/90')
    os.makedirs('../Resultados_similaridade/80')
    os.makedirs('../Resultados_similaridade/70')
    os.makedirs('../Resultados_similaridade/60')
    os.makedirs('../Codigos_busca_genes')
    os.makedirs('../Codigos_busca_genes/Codigos_genes_queries')
    os.makedirs('../Codigos_busca_genes/Listas_genes_pgp')
    os.makedirs('../Codigos_busca_genes/Listas_genes_pgp_faltantes')

# Função para realizar o Download dos genomas em formato gb
def download_genoma(genoma_id):
    
    saida = open('../Banco_de_dados/banco_de_dados_genomas/genoma_'+str(genoma_id)+'.gb','w')
    handle = Entrez.efetch(db='nuccore', id=genoma_id, rettype='gb')
    seqRecord = SeqIO.read(handle, format='gb')
    handle.close()
    saida.write(seqRecord.format('gb'))
    saida.close()

# Função para organizar os arquivos gb em fasta
def organizar_dados(genes_genoma_id):
    
    arq = open('../Banco_de_dados/banco_de_dados_genes/banco_de_dados_geral/dados_genes_'+str(genes_genoma_id)+'.fasta', 'w')
    arq_desconhecidos = open('../Banco_de_dados/banco_de_dados_genes/banco_de_dados_geral/dados_genes_desconhecidos_'+str(genes_genoma_id)+'.fasta', 'w')

    for rec in SeqIO.parse('../Banco_de_dados/banco_de_dados_genomas/genoma_'+str(genes_genoma_id)+'.gb', 'genbank'):
        
        org = ''
        
        for feature in rec.features:

            gene = ""
            header = ""
            aux = 0
            
            if (feature.type == 'source'):

                if ('strain' in feature.qualifiers):

                    org = feature.qualifiers['strain'].pop()

            if (feature.type == "CDS"):

                if ('translation' in feature.qualifiers):
                    translation = feature.qualifiers["translation"].pop()
                else:
                    translation = 'AUSENTE'
                
                if ('gene' in feature.qualifiers):

                    gene = feature.qualifiers['gene'].pop()

                else:
                    
                    inference = feature.qualifiers['inference'].pop()
                    indice = inference.find('RefSeq:')
                
                    if not indice == -1:
                        
                        reference_sequence = inference[(indice+7):]
                        handle = Entrez.efetch(db="protein", id=reference_sequence, rettype="gb")

                        for record in SeqIO.parse(handle,'gb'):

                            for feats in record.features:

                                if (feats.type == "CDS") and ("gene" in feats.qualifiers):

                                    gene += feats.qualifiers['gene'].pop() + ' '    
                        
                        if gene == "":

                            aux = 1
                            header = 'inference: ' + str(inference) + ' '

                    else:
                        
                        aux = 1
                        header = 'inference: ' + str(inference) + ' '

                for quali in feature.qualifiers.keys():

                    if quali != 'translation' and feature.qualifiers[quali]:

                        header += str(quali) + ": " + str(feature.qualifiers[quali].pop()) + " "

                if aux == 0:
                    
                    arq.write('>'+gene+' '+org+' '+header+'\n')
                    arq.write(translation+'\n')

                else:

                    arq_desconhecidos.write('>'+org+' '+header+'\n')
                    arq_desconhecidos.write(translation+'\n')

    arq.close()
    arq_desconhecidos.close()

# Função para criar o arquivo fasta dos genes encontrados anteriormente
def buscar_genes(genes_genoma_id):
    
    lista_genes = []
    arquivo = open('../Banco_de_dados/banco_de_dados_genes/banco_de_dados_pgp/genes_pgp_'+str(genes_genoma_id)+'.fasta','w')
    arquivo_final = open('../Banco_de_dados/banco_de_dados_blast/banco_de_dados_sbjcts/dados_genes_sbjct_'+str(genes_genoma_id)+'.fasta','a')

    with open('../Banco_de_dados/banco_de_dados_genes/banco_de_dados_geral/dados_genes_desconhecidos_'+str(genes_genoma_id)+'.fasta', 'r') as arquivo_sem_nome:

        arquivo_final.write(arquivo_sem_nome.read())


    for rec in SeqIO.parse('../Banco_de_dados/banco_de_dados_genes/banco_de_dados_geral/dados_genes_'+str(genes_genoma_id)+'.fasta', 'fasta'):

        lista_genes.append([str(rec.name),str(rec.description),str(rec.seq)])


    arq_genes = open('../Codigos_busca_genes/Listas_genes_pgp/genes_'+str(genes_genoma_id)+'.txt','r')
    lines = arq_genes.readlines()

    genes = []

    for line in lines:

        genes.append(line.strip())

    for gene_vet in lista_genes:

        gene = gene_vet[0]

        if gene in genes:

            arquivo.write('>'+gene_vet[1]+'\n')
            arquivo.write(gene_vet[2]+'\n')
            genes.remove(gene)

        else:

            arquivo_final.write('>'+gene_vet[1]+'\n')
            arquivo_final.write(gene_vet[2]+'\n')

    arquivo_faltantes = open('../Codigos_busca_genes/Listas_genes_pgp_faltantes/genes_pgp_faltantes_'+str(genes_genoma_id)+'.txt','w')
    for gen in genes:
        arquivo_faltantes.write(gen+'\n')
    arquivo_faltantes.close()
    arq_genes.close()
    arquivo.close()
    arquivo_final.close()
    arquivo_sem_nome.close()

# Leitura do Input do Script
if len(sys.argv) != 2:

    print("Uso: python download_e_organizacao_genomas.py input.txt")

else:

    nome_arquivo = sys.argv[1]
    valores = ler_arquivo(nome_arquivo)

    if valores:

        criar_pastas()

        for genoma_id in valores:

            download_genoma(genoma_id)
            organizar_dados(genes_genoma_id)
            buscar_genes(genes_genoma_id)

#Fim do Script
