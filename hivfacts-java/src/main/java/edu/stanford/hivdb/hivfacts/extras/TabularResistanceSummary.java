/*

    Copyright (C) 2017 Stanford HIVDB team

    Sierra is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Sierra is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

package edu.stanford.hivdb.hivfacts.extras;

import java.util.ArrayList;
import java.util.Collections;
import java.util.TreeMap;
import java.util.stream.Collectors;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;

import edu.stanford.hivdb.drugresistance.GeneDR;
import edu.stanford.hivdb.drugs.Drug;
import edu.stanford.hivdb.drugs.DrugClass;
import edu.stanford.hivdb.drugs.DrugResistanceAlgorithm;
import edu.stanford.hivdb.hivfacts.HIV;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Strain;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.mutations.MutationType;
import edu.stanford.hivdb.sequences.AlignedSequence;
import edu.stanford.hivdb.utilities.TSV;


/**
 * Print the sequence name, list of genes, "major" and "minor" mutations for each gene,
 *   resistance scores and levels for selected drugs
 * Currrently the header and list of drugs for a class (not the complete list) are constants
 *
 */
public class TabularResistanceSummary {
	private static final List<String> headerFields;

	static {
		HIV hiv = HIV.getInstance();

		List<String> hFields = new ArrayList<>();
		hFields.add("Sequence Name");
		hFields.add("Strain");
		hFields.add("Genes");
		for (DrugClass<HIV> drugClass : hiv.getDrugClasses()) {
			for (MutationType<HIV> mtype : drugClass.getMutationTypes()) {
				String mtypeText = mtype.getName();
				if (mtypeText.equals("Other")) {
					continue;
				}
				hFields.add(mtype.getFullName(drugClass));
			}
			for (Drug<HIV> drug : drugClass.getDrugs()) {
				hFields.add(drug.getDisplayAbbr() + " Score");
				hFields.add(drug.getDisplayAbbr() + " Level");
			}
		}
		hFields.add("Algorithm Name");
		hFields.add("Algorithm Version");
		hFields.add("Algorithm Date");
		headerFields = Collections.unmodifiableList(hFields);
	}

	private List<List<String>> sequenceRows = new ArrayList<>();
	private Map<String, Map<String, String>> tabularResults = new TreeMap<>();

	public TabularResistanceSummary(List<AlignedSequence<HIV>> alignedSeqs, List<Map<Gene<HIV>, GeneDR<HIV>>> allResistanceResults) {
		HIV hiv = HIV.getInstance();
		int numSeqs = alignedSeqs.size();
		for (int i=0; i < numSeqs; i ++) {
			AlignedSequence<HIV> alignedSeq = alignedSeqs.get(i);
			Map<Gene<HIV>, GeneDR<HIV>> resistanceResults = allResistanceResults.get(i);

			String seqName = alignedSeq.getInputSequence().getHeader();
			tabularResults.put(seqName, new TreeMap<String, String>());

			List<Gene<HIV>> geneList = Lists.newArrayList(resistanceResults.keySet());
			geneList.sort(Gene::compareTo);
			String genes = (
				geneList.stream()
				.map(gene -> gene.getAbstractGene())
				.collect(Collectors.joining(","))
			);
			Strain<HIV> strain = alignedSeq.getStrain();

			List<String> sequenceRecord = new ArrayList<>();
			sequenceRecord.add(seqName);
			sequenceRecord.add(strain.getName());
			sequenceRecord.add(genes);
			for (String absGene : hiv.getAbstractGenes()) {
				GeneDR<HIV> result = resistanceResults.getOrDefault(strain.getGene(absGene), null);
				sequenceRecord.addAll(getScoredMutations(absGene, result));
				sequenceRecord.addAll(getScores(absGene, result));
			}

			DrugResistanceAlgorithm<HIV> latest = hiv.getLatestDrugResistAlgorithm("HIVDB");
			sequenceRecord.add(latest.getFamily());
			sequenceRecord.add(latest.getVersion());
			sequenceRecord.add(latest.getPublishDate());
			sequenceRows.add(sequenceRecord);

			for (int j=0; j<headerFields.size(); j++) {
				String field = headerFields.get(j);
				String dataItem = sequenceRecord.get(j);
				tabularResults.get(seqName).put(field, dataItem);
			}
		}
	}

	@Override
	public String toString() {
		return TSV.dumps(headerFields, sequenceRows);
	}

	public Map<String, Map<String, String>> getTable() { return tabularResults; }
	public List<String> getHeaderFields() { return headerFields; }

	private static List<String> getScores(String absGene, GeneDR<HIV> geneDR) {
		HIV hiv = HIV.getInstance();
		List<String> resistanceScoresAndLevels = new ArrayList<>();
		for (DrugClass<HIV> drugClass : hiv.getDrugClasses()) {
			if (!drugClass.getAbstractGene().equals(absGene)) {
				continue;
			}
			for (Drug<HIV> drug : drugClass.getDrugs()) {
				if (geneDR == null) {
					resistanceScoresAndLevels.add("NA");
					resistanceScoresAndLevels.add("NA");
				}
				else {
					int score = geneDR.getTotalDrugScore(drug).intValue();
					int level = geneDR.getDrugLevel(drug);
					resistanceScoresAndLevels.add(Integer.toString(score));
					resistanceScoresAndLevels.add(Integer.toString(level));
				}
			}
		}
		return resistanceScoresAndLevels;
	}


	private static List<String> getScoredMutations(String absGene, GeneDR<HIV> geneDR) {
		HIV hiv = HIV.getInstance();
		List<String> scoredMutations = new ArrayList<>();
		Map<MutationType<HIV>, MutationSet<HIV>> mutationsByTypes;
		if (geneDR == null) {
			mutationsByTypes = Collections.emptyMap();
		}
		else {
			mutationsByTypes = geneDR.groupMutationsByTypes();
		}
		for (DrugClass<HIV> drugClass : hiv.getDrugClasses()) {
			if (!drugClass.getAbstractGene().equals(absGene)) {
				continue;
			}
			for (MutationType<HIV> mtype : drugClass.getMutationTypes()) {
				String mutationsText = "NA";
				if (mutationsByTypes.containsKey(mtype)) {
					MutationSet<HIV> muts = mutationsByTypes.get(mtype);
					if (muts.size() > 0) {
						mutationsText = muts.join();
					}
					else {
						mutationsText = "None";
					}
				}
				scoredMutations.add(mutationsText);
			}
		}
		return scoredMutations;
	}


}
