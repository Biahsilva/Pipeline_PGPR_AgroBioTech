from Bio import SeqIO
from Bio.Blast.Applications import*
from Bio.Blast import NCBIXML
import os
import sys

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

def blastp(genoma_id):

    arquivo_xml = open("../Relatorios_blast/relatorios_blast_complementar/relatorio_blastp_"+str(genoma_id)+".xml","r")
    dados = NCBIXML.parse(arquivo_xml)
    aux = 1
    arq = open('../Codigos_busca_genes/Listas_genes_pgp_faltantes/genes_pgp_faltantes_'+str(genoma_id)+'.txt')
    leitura = arq.readlines()
    arq.close()
    j=0
    rec = SeqIO.parse('../Banco_de_dados/banco_de_dados_blast/banco_de_dados_queries/query_'+str(genoma_id)+'.fasta','fasta')
    
    for item,sbj in zip(dados,rec):

        record = leitura[j].strip()
        j+=1
        i = 1

        if os.path.exists("../Relatorios_blast/relatorios_blast_complementar/relatorio_blastp_"+str(genoma_id)+".txt") and aux==1:
            os.remove("../Relatorios_blast/relatorios_blast_complementar/relatorio_blastp_"+str(genoma_id)+".txt")
        
        arquivo_final = open("../Relatorios_blast/relatorios_blast_complementar/relatorio_blastp_"+str(genoma_id)+".txt",'a')
        arquivo_final.write('ALINHAMENTO BLAST '+str(item.query)+'\n\n')

        if os.path.exists("../Banco_de_dados/banco_de_dados_genes/banco_de_dados_pgp/genes_pgp_complementares_"+str(genoma_id)+".fasta") and aux==1:
            os.remove("../Banco_de_dados/banco_de_dados_genes/banco_de_dados_pgp/genes_pgp_complementares_"+str(genoma_id)+".fasta")
        
        arquivo_complementar = open("../Banco_de_dados/banco_de_dados_genes/banco_de_dados_pgp/genes_pgp_complementares_"+str(genoma_id)+".fasta","a")

        for a in item.alignments:

            aux = 0

            for hsp in a.hsps:

                if hsp.identities/a.length > 0.9:

                    arquivo_final.write('Alinhamento '+str(i)+'\n')
                    arquivo_final.write('Sequencia: '+str(a.title)+'\n')
                    arquivo_final.write('Tamanho: '+str(a.length)+'\n')
                    arquivo_final.write('Score: '+str(hsp.score)+'\n')
                    arquivo_final.write('Bits: '+str(hsp.bits)+'\n')
                    percentage = "{:.0%}".format(hsp.identities/a.length)
                    arquivo_final.write('Identities: '+str(percentage)+'\n')
                    arquivo_final.write('E-Value: '+str(hsp.expect)+'\n')
                    arquivo_final.write(str(hsp.query)+'\n')
                    arquivo_final.write(str(hsp.match)+'\n')
                    arquivo_final.write(str(hsp.sbjct)+'\n')
                    arquivo_final.write("\n\n\n")
                    arquivo_complementar.write('>'+str(record)+' '+str(a.title)+'\n')
                    arquivo_complementar.write(str(sbj.seq)+'\n')

                else:

                    arquivo_final.write('\n\nO ALINHAMENTO NÃO DEU MAIS DE 90% DE IDENTIDADE!!\n\n')
                
                i+=1
        
        arquivo_final.write('------------------------------------------------------------------------------------------------\n\n\n')

# Leitura do Input do Script
if len(sys.argv) != 2:

    print("Uso: python analisar_blastp.py input.txt")

else:

    nome_arquivo = sys.argv[1]
    valores = ler_arquivo(nome_arquivo)

    if valores:

        for genoma_id in valores:
            
            blastp(genoma_id)

#Fim do Script
