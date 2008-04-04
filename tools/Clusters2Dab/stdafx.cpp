#include "stdafx.h"

/*!
 * \page Clusters2Dab Clusters2Dab
 * 
 * Clusters2Dab converts the output of many common non-hierarchical clustering algorithms into pairwise
 * scores based on the frequency or confidence of coclustering.  For example, if two genes are scored by
 * the number of times they cocluster over many random seeds, then high scores will be indicative of a
 * stronger pairwise relationship than low scores.
 * 
 * \section sec_overview Overview
 * 
 * \subsection ssec_usage_detailed Detailed Usage
 * 
 * \include Clusters2Dab/Clusters2Dab.ggo
 * 
 * <table><tr>
 *	<th>Flag</th>
 *	<th>Default</th>
 *	<th>Type</th>
 *	<th>Description</th>
 * </tr><tr>
 *	<td>-i</td>
 *	<td>stdin</td>
 *	<td>Text file</td>
 *	<td>Input cluster file in one of the text formats supported by \c -t.</td>
 * </tr><tr>
 *	<td>-o</td>
 *	<td>stdout</td>
 *	<td>DAT/DAB file</td>
 *	<td>Output DAT/DAB file containing pairwise scores appropriate to the input cluster type.</td>
 * </tr><tr>
 *	<td>-t</td>
 *	<td>list</td>
 *	<td>list, samba, param, or fuzzy</td>
 *	<td>Type of cluster file provided to \c -i.  \c list assumes that each line contains a gene ID and a
 *		cluster index separated by a tab, \c samba reads biclustering output from the EXPANDER program by
 *		Sharan, Shamir, et al, \c param reads hard clustering output from EXPANDER, and \c fuzzy reads
 *		output from the Aerie fuzzy k-means program by Gasch, Eisen, et al.</td>
 * </tr><tr>
 *	<td>-c</td>
 *	<td>off</td>
 *	<td>Flag</td>
 *	<td>If on, calculate pairwise scores solely by cocluster frequency (counts); otherwise, pairwise scores
 *		are weighted by cluster confidence.</td>
 * </tr><tr>
 *	<td>-s</td>
 *	<td>2</td>
 *	<td>Integer</td>
 *	<td>Number of columns to skip between the initial ID column and the first experimental (data) column
 *		in the input PCL.</td>
 * </tr></table>
 */
