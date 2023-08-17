from Bio import SeqIO


def concatenar_genes(genoma_id):

    arquivo_pgp = open('../Banco_de_dados/banco_de_dados_genes/banco_de_dados_pgp/genes_pgp_'+str(genoma_id)+'.fasta','r')
    arquivo_pgp_complementar = open('../Banco_de_dados/banco_de_dados_genes/banco_de_dados_pgp/genes_pgp_complementares_'+str(genoma_id)+'.fasta','r')
    arquivo_pgp_AgroBioTech = open('../Banco_de_dados/banco_de_dados_AgroBioTech/genes_pgp_referencias/genes_pgp_referencias_'+str(genoma_id)+'.fasta','w')

    arquivo_pgp_AgroBioTech.write(arquivo_pgp.read())
    arquivo_pgp_AgroBioTech.write(arquivo_pgp_complementar.read())

    arquivo_pgp_AgroBioTech.close()
    arquivo_pgp_complementar.close()
    arquivo_pgp.close()

def minerar_genes(genoma_id):
    
    genes = []

    for rec in SeqIO.parse('../Banco_de_dados/banco_de_dados_AgroBioTech/genes_pgp_referencias/genes_pgp_referencias_'+str(genoma_id)+'.fasta', "fasta"):
        
        genes.append(rec.name)

    return genes

lista_id_ref = ['CP040881.1', 'CP029034.1', 'CP099465.1', 'CP000560.2']
lista_genomas_ref = ['HNA3', 'LDO2', 'AG75','FZB42']

lista_id_labim = ['CP045993.1', 'CP023748.1', 'CP079719.1']
lista_genomas_labim = ['LABIM22', 'LABIM40', 'LABIM44']

for genoma_id in lista_id_ref+lista_id_labim:
    
    concatenar_genes(genoma_id)


dicionario_genomas_ref = {}
dicionario_genomas_labim = {}

for genoma_id_ref,genoma_name_ref in zip(lista_id_ref,lista_genomas_ref):
    dicionario_genomas_ref[genoma_name_ref] = set(minerar_genes(genoma_id_ref))
    
for genoma_id_labim,genoma_name_labim in zip(lista_id_labim,lista_genomas_labim):
    dicionario_genomas_labim[genoma_name_labim] = set(minerar_genes(genoma_id_labim))


# HNA3 - CP040881.1
# LDO2 - CP029034.1
# AG75 - CP099465.1
# RZ2MS9 - CP049978.1 #Única que não usa
# FZB42 - CP000560.2
# LABIM22 - CP045993.1
# LABIM40 - CP023748.1
# LABIM44 - CP079719.1

teste=set([])

for i in lista_genomas_ref:

    teste.update(dicionario_genomas_ref[i])

print('tamanho '+str(len(teste))+'\n\n')




arquivo_relatorio = open('../Resultados_comparacao/resultado_comparacao.txt','w')

arquivo_relatorio.write('../Relatório de Comparação dos Genes do Banco de Dados AgroBioTech - Parte 1\n\n')

#Interseção 1 - LABIM's contra Referências

for chaves_labim in dicionario_genomas_labim.keys():

    for chaves_ref in dicionario_genomas_ref.keys():

        print('Comparação entre os Genomas '+chaves_labim+' e '+chaves_ref+'.\n')
        intersect = dicionario_genomas_labim[chaves_labim].intersection(dicionario_genomas_ref[chaves_ref])
        arquivo_relatorio.write('Comparação entre os Genomas '+chaves_labim+' e '+chaves_ref+': '+str(len(intersect))+'\n'+str(intersect)+'.\n\n')


arquivo_relatorio.write('\n\nRelatório de Comparação dos Genes do Banco de Dados AgroBioTech - Parte 2\n\n')



#Interseção 2 - Referências entre si

for chaves_ref1 in dicionario_genomas_ref.keys():

    for chaves_ref2 in dicionario_genomas_ref.keys():

        if not chaves_ref1==chaves_ref2:

            print('Comparação entre os Genomas '+chaves_ref1+' e '+chaves_ref2+'.\n')
            intersect = dicionario_genomas_ref[chaves_ref1].intersection(dicionario_genomas_ref[chaves_ref2])
            arquivo_relatorio.write('Comparação entre os Genomas '+chaves_ref1+' e '+chaves_ref2+': '+str(len(intersect))+'\n'+str(intersect)+'.\n\n')

arquivo_relatorio.write('\n\n')


for chaves_ref1 in dicionario_genomas_ref.keys():

    for chaves_ref2 in dicionario_genomas_ref.keys():

        for chaves_ref3 in dicionario_genomas_ref.keys():
            
            if not (chaves_ref1==chaves_ref2 or chaves_ref1==chaves_ref3 or chaves_ref2==chaves_ref3):

                print('Comparação entre os Genomas '+chaves_ref1+', '+chaves_ref2+' e '+chaves_ref3+'.\n')
                intersect = dicionario_genomas_ref[chaves_ref1].intersection(dicionario_genomas_ref[chaves_ref2].intersection(dicionario_genomas_ref[chaves_ref3]))
                arquivo_relatorio.write('Comparação entre os Genomas '+chaves_ref1+', '+chaves_ref2+' e '+chaves_ref3+': '+str(len(intersect))+'\n'+str(intersect)+'.\n\n')

intersect = dicionario_genomas_ref['HNA3'].intersection(dicionario_genomas_ref['AG75'].intersection(dicionario_genomas_ref['LDO2'].intersection(dicionario_genomas_ref['FZB42'])))
arquivo_relatorio.write('Comparação entre os Genomas HNA3, AG75, LDO2 e FZB42: '+str(len(intersect))+'\n'+str(intersect)+'.\n\n')

arquivo_relatorio.close()
