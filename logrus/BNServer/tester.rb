#!/usr/bin/ruby

require "socket"

sock = TCPSocket.new( "localhost", 1234 )

=begin # Perform inference
# size, opcode, context, genes
sock.write( [9, 0, 1, 0].pack( "ICII" ) )
iBytes = sock.read( 4 ).unpack( "I" )[ 0 ]
sock.read( iBytes ).unpack( "f" * ( iBytes / 4 ) ).each do |d|
	puts( d ); end
=end

=begin # Request data
# size, opcode, gene1, gene2
sock.write( [9, 1, 1, 4].pack( "ICII" ) )
iBytes = sock.read( 4 ).unpack( "I" )[ 0 ]
sock.read( iBytes ).unpack( "C" * iBytes ).each do |i|
	puts( i ); end
=end

#=begin # Pixie graph
# size, opcode, context, limit, genes
sock.write( [13, 2, 0, 3, 1].pack( "ICIII" ) )
iBytes = sock.read( 4 ).unpack( "I" )[ 0 ]
sock.read( iBytes ).unpack( "C" * iBytes ).each do |i|
	puts( i ); end
#=end
