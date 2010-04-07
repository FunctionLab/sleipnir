#!/usr/bin/ruby

if( ARGV.length != 1 )
	raise( "Usage: half2weights.rb <bn.xdsl> < <mi_half.txt>" ); end
strBN = ARGV[ 0 ]

astrIDs = []
IO.foreach( strBN ) do |strLine|
	if( ( strLine =~ /cpt id=\"([^"]+)\"/ ) && ( $1 != "FR" ) ) # "
		astrIDs.push( $1 ); end; end

astrHeaders = adSums = nil
adSelf = []
STDIN.each do |strLine|
	astrLine = strLine.strip.split( /\t/ )
	if( !astrHeaders )
		astrHeaders = astrLine.map do |str|
			str =~ /^(?:.*\/)?(\S+)\.dab$/
			$1; end
		adSums = Array.new( astrHeaders.length, 0 )
		next; end
	adSelf[ $. - 2 ] = astrLine[ $. - 1 ].to_f
	($.).upto( astrLine.length - 1 ) do |i|
		d = astrLine[ i ].to_f
		adSums[ $. - 2 ] += d
		adSums[ i - 1 ] += d; end; end

hashWeights = {}
astrHeaders.each_index do |i|
	hashWeights[ astrHeaders[ i ] ] =
# Math.log( adSums[ i ] / adSelf[ i ] + 1 ) # log
# adSums[ i ] / adSelf[ i ]
2 ** ( adSums[ i ] / adSelf[ i ] ) - 1
end

astrIDs.each do |strID|
	puts( [strID, hashWeights[ strID ] || "0"].join( "\t" ) ); end
