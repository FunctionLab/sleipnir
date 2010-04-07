#!/usr/bin/ruby

astrHeaders = nil
aadData = []
STDIN.each do |strLine|
	astrLine = strLine.split( /\t/ ).map do |str|
		str.strip; end
	if( !astrHeaders )
		astrHeaders = astrLine
		next; end
	1.upto( $. - 2 ) do |i|
		astrLine[ i ] = aadData[ i - 1 ][ $. - 2 ]; end
	aadData.push( astrLine[ 1, astrLine.length ].map do |str|
		str.to_f; end ); end

aadData.each_index do |i|
	aadData[ i ].each_index do |j|
		if( i == j )
			next; end
		if( !( aadData[ i ][ i ] && aadData[ j ][ j ] ) )
			aadData[ i ][ j ] = nil
		else
			if( ( d = ( aadData[ i ][ j ] / [aadData[ i ][ i ],
				aadData[ j ][ j ]].min ) ) > 1 )
				raise( "Invalid value " + d.to_s + " >1: " +
					aadData[ i ][ j ].to_s + " / " + [aadData[ i ][ i ],
					aadData[ j ][ j ]].inspect + ".min for " +
					astrHeaders[ i + 1 ] + ", " + astrHeaders[ j + 1 ] ); end
			aadData[ i ][ j ] = d; end; end; end
aadData.each_index do |i|
	aadData[ i ][ i ] = 1; end

puts( astrHeaders.join( "\t" ) )
aadData.each_index do |i|
	puts( aadData[ i ].unshift( astrHeaders[ i + 1 ] ).join( "\t" ) ); end
