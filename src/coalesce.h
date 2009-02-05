/*****************************************************************************
* This file is provided under the Creative Commons Attribution 3.0 license.
*
* You are free to share, copy, distribute, transmit, or adapt this work
* PROVIDED THAT you attribute the work to the authors listed below.
* For more information, please see the following web page:
* http://creativecommons.org/licenses/by/3.0/
*
* This file is a component of the Sleipnir library for functional genomics,
* authored by:
* Curtis Huttenhower (chuttenh@princeton.edu)
* Mark Schroeder
* Maria D. Chikina
* Olga G. Troyanskaya (ogt@princeton.edu, primary contact)
*
* If you use this library, the included executable tools, or any related
* code in your work, please cite the following publication:
* Curtis Huttenhower, Mark Schroeder, Maria D. Chikina, and
* Olga G. Troyanskaya.
* "The Sleipnir library for computational functional genomics"
*****************************************************************************/
#ifndef COALESCE_H
#define COALESCE_H

#include "coalescei.h"
#include "coalescecluster.h"
#include "coalescemotifs.h"

namespace Sleipnir {

class CCoalesce : CCoalesceImpl {
public:
	bool Cluster( const CPCL& PCL, const CFASTA& FASTA, std::vector<CCoalesceCluster>& vecClusters );

	void SetPValueCorrelation( float dPValue ) {

		m_dPValueCorrelation = dPValue; }

	float GetPValueCorrelation( ) const {

		return m_dPValueCorrelation; }

	void SetBins( size_t iBins ) {

		m_iBins = iBins; }

	size_t GetBins( ) const {

		return m_iBins; }

	float GetPValueCondition( ) const {

		return m_dPValueCondition; }

	void SetPValueCondition( float dPValue ) {

		m_dPValueCondition = dPValue; }

	float GetPValueMotif( ) const {

		return m_dPValueMotif; }

	void SetPValueMotif( float dPValue ) {

		m_dPValueMotif = dPValue; }

	float GetProbabilityGene( ) const {

		return m_dProbabilityGene; }

	void SetProbabilityGene( float dProbability ) {

		m_dProbabilityGene = dProbability; }

	bool IsOutputIntermediate( ) const {

		return !GetDirectoryIntermediate( ).empty( ); }

	const std::string& GetDirectoryIntermediate( ) const {

		return m_strDirectoryIntermediate; }

	void SetDirectoryIntermediate( const std::string& strDirectoryIntermediate ) {

		m_strDirectoryIntermediate = strDirectoryIntermediate; }

	void SetMotifs( CCoalesceMotifLibrary& Motifs ) {

		if( m_fMotifs && m_pMotifs && ( m_pMotifs != &Motifs ) )
			delete m_pMotifs;
		m_pMotifs = &Motifs; }

	const CCoalesceMotifLibrary* GetMotifs( ) const {

		return m_pMotifs; }

	size_t GetK( ) const {

		return m_iK; }

	void SetK( size_t iK ) {

		m_iK = iK; }

	size_t GetBasesPerMatch( ) const {

		return m_iBasesPerMatch; }

	void SetBasesPerMatch( size_t iBasesPerMatch ) {

		m_iBasesPerMatch = iBasesPerMatch; }

	const std::string& GetSequenceCache( ) const {

		return m_strSequenceCache; }

	void SetSequenceCache( const std::string& strSequenceCache ) {

		m_strSequenceCache = strSequenceCache; }

	float GetPValueMerge( ) const {

		return m_dPValueMerge; }

	void SetPValueMerge( float dPValue ) {

		m_dPValueMerge = dPValue; }

	float GetCutoffMerge( ) const {

		return m_dCutoffMerge; }

	void SetCutoffMerge( float dCutoff ) {

		m_dCutoffMerge = dCutoff; }

	size_t GetSizeMinimum( ) const {

		return m_iSizeMinimum; }

	void SetSizeMinimum( size_t iSizeGenes ) {

		m_iSizeMinimum = iSizeGenes; }

	size_t GetSizeMaximum( ) const {

		return m_iSizeMaximum; }

	void SetSizeMaximum( size_t iSizeMotifs ) {

		m_iSizeMaximum = iSizeMotifs; }

	size_t GetSizeMerge( ) const {

		return m_iSizeMerge; }

	void SetSizeMerge( size_t iSizeMerge ) {

		m_iSizeMerge = iSizeMerge; }

	void ClearDatasets( ) {

		m_vecsDatasets.clear( ); }

	bool AddDataset( const std::set<size_t>& setiDataset ) {
		size_t								i;
		std::set<size_t>::const_iterator	iterExperiment;

		if( setiDataset.empty( ) )
			return true;
		for( iterExperiment = setiDataset.begin( ); iterExperiment != setiDataset.end( ); ++iterExperiment )
			for( i = 0; i < m_vecsDatasets.size( ); ++i )
				if( m_vecsDatasets[ i ].IsCondition( *iterExperiment ) )
					return false;

		m_vecsDatasets.push_back( SCoalesceDataset( setiDataset ) );
		return true; }

	void SetNumberCorrelation( size_t iPairs ) {

		m_iNumberCorrelation = iPairs; }

	size_t GetNumberCorrelation( ) const {

		return m_iNumberCorrelation; }

	void SetThreads( size_t iThreads ) {

		m_iThreads = iThreads; }

	size_t GetThreads( ) const {

		return m_iThreads; }

	void AddWiggle( const CFASTA& FASTA ) {

		m_vecpWiggles.push_back( &FASTA ); }

	void ClearWiggles( ) {

		m_vecpWiggles.clear( ); }
};

}

#endif // COALESCE_H
