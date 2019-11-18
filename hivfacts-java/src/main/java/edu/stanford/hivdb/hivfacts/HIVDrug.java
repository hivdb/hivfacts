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

package edu.stanford.hivdb.hivfacts;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

import edu.stanford.hivdb.hivfacts.HIVDrugClass;

public enum HIVDrug implements Drug<HIVDrugClass> {
	ABC("abacavir", "ABC", HIVDrugClass.NRTI),
	AZT("zidovudine", "AZT", HIVDrugClass.NRTI),
	D4T("stavudine", "D4T",HIVDrugClass.NRTI),
	DDI("didanosine", "DDI", HIVDrugClass.NRTI),
	FTC("emtricitabine", "FTC", HIVDrugClass.NRTI),
	LMV("lamivudine", "3TC", HIVDrugClass.NRTI),
	TDF("tenofovir", "TDF", HIVDrugClass.NRTI),

	ATV("atazanavir/r", "ATV/r", HIVDrugClass.PI),
	DRV("darunavir/r", "DRV/r", HIVDrugClass.PI),
	FPV("fosamprenavir/r", "FPV/r", HIVDrugClass.PI),
	IDV("indinavir/r", "IDV/r", HIVDrugClass.PI),
	LPV("lopinavir/r", "LPV/r", HIVDrugClass.PI),
	NFV("nelfinavir", "NFV", HIVDrugClass.PI),
	SQV("saquinavir/r", "SQV/r", HIVDrugClass.PI),
	TPV("tipranavir/r", "TPV/r", HIVDrugClass.PI),

	DOR("doravirine", "DOR", HIVDrugClass.NNRTI),
	EFV("efavirenz", "EFV", HIVDrugClass.NNRTI),
	ETR("etravirine", "ETR", HIVDrugClass.NNRTI),
	NVP("nevirapine", "NVP", HIVDrugClass.NNRTI),
	RPV("rilpivirine", "RPV", HIVDrugClass.NNRTI),

	BIC("bictegravir", "BIC", HIVDrugClass.INSTI),
	DTG("dolutegravir", "DTG", HIVDrugClass.INSTI),
	EVG("elvitegravir", "EVG", HIVDrugClass.INSTI),
	RAL("raltegravir", "RAL",HIVDrugClass.INSTI);

	private static final Map<String, HIVDrug> drugSynonyms;

	private final String fullName;
	private final String displayAbbr;
	private final HIVDrugClass drugClass;

	static {
		Map<String, HIVDrug> tmpDrugSynonyms = new HashMap<>();
		tmpDrugSynonyms.put("LPV/r", HIVDrug.LPV);
		tmpDrugSynonyms.put("IDV/r", HIVDrug.IDV);
		tmpDrugSynonyms.put("FPV/r", HIVDrug.FPV);
		tmpDrugSynonyms.put("ATV/r", HIVDrug.ATV);
		tmpDrugSynonyms.put("DRV/r", HIVDrug.DRV);
		tmpDrugSynonyms.put("TPV/r", HIVDrug.TPV);
		tmpDrugSynonyms.put("SQV/r", HIVDrug.SQV);
		tmpDrugSynonyms.put("3TC",  HIVDrug.LMV);
		// ANRS distinguishes QD and BID for DRV/r and DTG;
		// By default we only use BID here
		tmpDrugSynonyms.put("DRV/r_QD", HIVDrug.DRV);
		tmpDrugSynonyms.put("DTG_QD", HIVDrug.DTG);
		drugSynonyms = Collections.unmodifiableMap(tmpDrugSynonyms);
	}

	private HIVDrug(
		final String fullName, final String displayAbbr,
		final HIVDrugClass drugClass
	) {
		this.fullName = fullName;
		this.displayAbbr = displayAbbr;
		this.drugClass = drugClass;
	}

	@Override
	public HIVDrugClass getDrugClass() {
		return drugClass;
	}

	@Override
	public String getFullName() {
		return fullName;
	}

	@Override
	public String getDisplayAbbr() {
		return displayAbbr;
	}

	public static HIVDrug getSynonym(String synonym) {
		if (drugSynonyms.containsKey(synonym)) {
			return drugSynonyms.get(synonym);
		} else {
			try {
				return HIVDrug.valueOf(synonym);
			} catch (IllegalArgumentException e) {
				return null;
			}
		}
	}
}
