#!/usr/local/bin/ruby

require "socket"

sock = TCPSocket.new( "localhost", 3009 )

=begin # Perform inference
# size, opcode, context, genes
sock.write( [9, 0, 11, 1].pack( "ICII" ) )
iBytes = sock.read( 4 ).unpack( "I" )[ 0 ]
sock.read( iBytes ).unpack( "f" * ( iBytes / 4 ) ).each do |d|
	puts( d ); end
=end

=begin # Request data
# size, opcode, gene1, gene2
sock.write( [9, 1, 1, 3].pack( "ICII" ) )
iBytes = sock.read( 4 ).unpack( "I" )[ 0 ]
sock.read( iBytes ).unpack( "C" * iBytes ).each do |i|
	puts( i ); end
=end

=begin # Pixie graph
# size, opcode, file, context, limit, filter, genes
sock.write( [18, 2, 0, 0, 3, 1, 1].pack( "ICCIIFI" ) )
iBytes = sock.read( 4 ).unpack( "I" )[ 0 ]
puts( sock.read( iBytes ) )
=end

=begin # Contexts
# size, opcode, genes
sock.write( [5, 3, 1].pack( "ICI" ) )
iBytes = sock.read( 4 ).unpack( "I" )[ 0 ]
sock.read( iBytes ).unpack( "f" * ( iBytes / 4 ) ).each do |d|
	puts( d ); end
=end

=begin # Diseases
# size, opcode, context, genes
sock.write( [9, 5, 0, 1].pack( "ICII" ) )
iBytes = sock.read( 4 ).unpack( "I" )[ 0 ]
sock.read( iBytes ).unpack( "f" * ( iBytes / 4 ) ).each do |d|
	puts( d ); end
=end

=begin # Genes
# size, opcode, context, genes
sock.write( [17, 6, 0, 1, 2, 3].pack( "ICIIII" ) )
iBytes = sock.read( 4 ).unpack( "I" )[ 0 ]
sock.read( iBytes ).unpack( "f" * ( iBytes / 4 ) ).each do |d|
	puts( d ); end
=end

=begin # Association
# size, opcode, context, genes1, 0, genes2
sock.write( [33, 7, 0, 1, 2, 3, 0, 4, 5, 6].pack( "ICIIIIIIII" ) )
iBytes = sock.read( 4 ).unpack( "I" )[ 0 ]
sock.read( iBytes ).unpack( "f" * ( iBytes / 4 ) ).each do |d|
	puts( d ); end
=end

=begin # Associations
# size, opcode, diseases?, context, genes
C_aiCarbohydrateTransport = [2456, 3435, 4578, 6289, 6418, 6567, 6918, 8474, 8690, 8818, 9396, 10091, 10233, 10280, 11280, 13457, 14776, 14910, 16406, 16409, 16418, 16427, 16437, 16447, 16459, 16471, 17927, 18150, 18445, 18459, 18469, 18471, 18485, 18656, 18750, 19526, 19534, 19536, 21093, 21100, 21387]
sock.write( [6 + ( 4 * C_aiCarbohydrateTransport.length ), 8, 0, -1].concat(
	C_aiCarbohydrateTransport ).pack( "ICCI" +
	( "I" * C_aiCarbohydrateTransport.length ) ) )
iBytes = sock.read( 4 ).unpack( "I" )[ 0 ]
sock.read( iBytes ).unpack( "f" * ( iBytes / 4 ) ).each do |d|
	puts( d ); end
=end
