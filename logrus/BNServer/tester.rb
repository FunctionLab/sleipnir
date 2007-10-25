#!/usr/bin/ruby

require "socket"

sock = TCPSocket.new( "localhost", 1234 )

=begin # Perform inference
sock.write( [9, 0, 1, 0].pack( "ICII" ) )
iBytes = sock.read( 4 ).unpack( "I" )[ 0 ]
sock.read( iBytes ).unpack( "f" * ( iBytes / 4 ) ).each do |d|
	puts( d ); end
=end

#=begin # Request data
sock.write( [9, 1, 1, 4].pack( "ICII" ) )
iBytes = sock.read( 4 ).unpack( "I" )[ 0 ]
sock.read( iBytes ).unpack( "C" * iBytes ).each do |i|
	puts( i ); end
#=end
