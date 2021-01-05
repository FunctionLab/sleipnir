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

genes = ['90634', '23659']
# genes = ['121212', '131313']
params = QueryParams()
query = SeekQuery(species='mock', genes=genes, parameters=params)
result = client.seek_query(query)
if result.success is True:
    for gs in result.gene_scores:
        print(f'gene: {gs.name}, {gs.value}')
    for ds in result.dataset_weights:
        print(f'dset: {ds.name}, {ds.value}')
else:
    print(f'query error: {result.statusMsg}')

transport.close()