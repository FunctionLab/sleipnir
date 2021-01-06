from thrift.transport import TTransport, TSocket
from thrift.protocol.TBinaryProtocol import TBinaryProtocol

from pyseek import SeekRPC
from pyseek.ttypes import SeekQuery, QueryParams, QueryResult

host = 'localhost'
port = 9090

socket = TSocket.TSocket(host, port)
transport = TTransport.TBufferedTransport(socket)
protocol = TBinaryProtocol(transport)
client = SeekRPC.Client(protocol)
transport.open()

#species = 'human'
#genes = ['90634', '23659']

species = 'fly'
genes = ['35234', '35232']
#genes = ['34930', '35234', '35232']

#species = 'yeast'
#genes = ['YGL142C', 'YHR188C']


params = QueryParams(distance_measure="ZscoreHubbinessCorrected",
                     min_query_genes_fraction=0.5,
                     min_genome_fraction=0.5)
query = SeekQuery(species=species, genes=genes, parameters=params)
result = client.seek_query(query)
if result.success is True:
    count = 0;
    for gs in result.gene_scores:
        print(f'gene: {gs.name}, {gs.value}')
        count += 1
        if count > 100: break

    count = 0;
    for ds in result.dataset_weights:
        print(f'dset: {ds.name}, {ds.value}')
        count += 1
        if count > 100: break
else:
    print(f'query error: {result.statusMsg}')

transport.close()
