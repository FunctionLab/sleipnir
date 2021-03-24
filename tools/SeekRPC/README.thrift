# Generate thrift interface files
thrift -gen cpp seek_rpc.thrift
thrift -gen py seek_rpc.thrift

cd gen-cpp
g++ --std=c++11 -c *.cpp
cd ..

# compile server
g++ --std=c++11 -Igen-cpp/ -o thriftserver thriftserver.cpp SeekInterface.cpp gen-cpp/SeekRPC.o gen-cpp/seek*.o -lthrift

# compile client
g++ --std=c++11 -Igen-cpp/ -o SeekRPCClient SeekRPCClient.cpp gen-cpp/seek*.o gen-cpp/SeekRPC.o -lthrift
g++ --std=c++11 -Igen-cpp/ -o thriftclient thriftclient.cpp gen-cpp/seek*.o gen-cpp/SeekRPC.o -lthrift

# old compile server
g++ --std=c++11 -Igen-cpp/ -o thriftserver gen-cpp/SeekRPC_server.skeleton.cpp gen-cpp/SeekRPC.o gen-cpp/seek*.o -lthrift
g++ --std=c++11 -Igen-cpp/ -o thriftserver server.cpp gen-cpp/SeekRPC.o gen-cpp/seek*.o -lthrift

# old compile client
g++ --std=c++11 -Igen-cpp/ -o thriftclient client.cpp gen-cpp/seek*.o gen-cpp/SeekRPC.o -lthrift

