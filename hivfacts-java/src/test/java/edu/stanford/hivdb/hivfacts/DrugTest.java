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

import static org.junit.Assert.*;

import org.junit.Test;

public class DrugTest {

	@Test
	public void testDrugCount() {
		// This test will fail every time we add a new drug.
		// WARNING: Don't just fix this test case. Add the new drug
		// to all following cases too.
		assertEquals(HIVDrug.values().length, 24);
	}

	@Test
	public void testGetDrugClass() {
		final HIVDrugClass PI = HIVDrugClass.PI;
		final HIVDrugClass NRTI = HIVDrugClass.NRTI;
		final HIVDrugClass NNRTI = HIVDrugClass.NNRTI;
		final HIVDrugClass INSTI = HIVDrugClass.INSTI;

		assertEquals(HIVDrug.ABC.getDrugClass(), NRTI);
		assertEquals(HIVDrug.AZT.getDrugClass(), NRTI);
		assertEquals(HIVDrug.D4T.getDrugClass(), NRTI);
		assertEquals(HIVDrug.DDI.getDrugClass(), NRTI);
		assertEquals(HIVDrug.FTC.getDrugClass(), NRTI);
		assertEquals(HIVDrug.LMV.getDrugClass(), NRTI);
		assertEquals(HIVDrug.TDF.getDrugClass(), NRTI);

		assertEquals(HIVDrug.ATV.getDrugClass(), PI);
		assertEquals(HIVDrug.DRV.getDrugClass(), PI);
		assertEquals(HIVDrug.FPV.getDrugClass(), PI);
		assertEquals(HIVDrug.IDV.getDrugClass(), PI);
		assertEquals(HIVDrug.LPV.getDrugClass(), PI);
		assertEquals(HIVDrug.NFV.getDrugClass(), PI);
		assertEquals(HIVDrug.SQV.getDrugClass(), PI);
		assertEquals(HIVDrug.TPV.getDrugClass(), PI);

		assertEquals(HIVDrug.EFV.getDrugClass(), NNRTI);
		assertEquals(HIVDrug.ETR.getDrugClass(), NNRTI);
		assertEquals(HIVDrug.NVP.getDrugClass(), NNRTI);
		assertEquals(HIVDrug.RPV.getDrugClass(), NNRTI);

		assertEquals(HIVDrug.DTG.getDrugClass(), INSTI);
		assertEquals(HIVDrug.EVG.getDrugClass(), INSTI);
		assertEquals(HIVDrug.RAL.getDrugClass(), INSTI);
	}

	@Test
	public void testGetFullName() {
		assertEquals(HIVDrug.ABC.getFullName(), "abacavir");
		assertEquals(HIVDrug.AZT.getFullName(), "zidovudine");
		assertEquals(HIVDrug.D4T.getFullName(), "stavudine");
		assertEquals(HIVDrug.DDI.getFullName(), "didanosine");
		assertEquals(HIVDrug.FTC.getFullName(), "emtricitabine");
		assertEquals(HIVDrug.LMV.getFullName(), "lamivudine");
		assertEquals(HIVDrug.TDF.getFullName(), "tenofovir");

		assertEquals(HIVDrug.ATV.getFullName(), "atazanavir/r");
		assertEquals(HIVDrug.DRV.getFullName(), "darunavir/r");
		assertEquals(HIVDrug.FPV.getFullName(), "fosamprenavir/r");
		assertEquals(HIVDrug.IDV.getFullName(), "indinavir/r");
		assertEquals(HIVDrug.LPV.getFullName(), "lopinavir/r");
		assertEquals(HIVDrug.NFV.getFullName(), "nelfinavir");
		assertEquals(HIVDrug.SQV.getFullName(), "saquinavir/r");
		assertEquals(HIVDrug.TPV.getFullName(), "tipranavir/r");

		assertEquals(HIVDrug.EFV.getFullName(), "efavirenz");
		assertEquals(HIVDrug.ETR.getFullName(), "etravirine");
		assertEquals(HIVDrug.NVP.getFullName(), "nevirapine");
		assertEquals(HIVDrug.RPV.getFullName(), "rilpivirine");

		assertEquals(HIVDrug.DTG.getFullName(), "dolutegravir");
		assertEquals(HIVDrug.EVG.getFullName(), "elvitegravir");
		assertEquals(HIVDrug.RAL.getFullName(), "raltegravir");
	}

	@Test
	public void testGetDisplayAbbr() {
		assertEquals(HIVDrug.ABC.getDisplayAbbr(), "ABC");
		assertEquals(HIVDrug.AZT.getDisplayAbbr(), "AZT");
		assertEquals(HIVDrug.D4T.getDisplayAbbr(), "D4T");
		assertEquals(HIVDrug.DDI.getDisplayAbbr(), "DDI");
		assertEquals(HIVDrug.FTC.getDisplayAbbr(), "FTC");
		assertEquals(HIVDrug.LMV.getDisplayAbbr(), "3TC");
		assertEquals(HIVDrug.TDF.getDisplayAbbr(), "TDF");

		assertEquals(HIVDrug.ATV.getDisplayAbbr(), "ATV/r");
		assertEquals(HIVDrug.DRV.getDisplayAbbr(), "DRV/r");
		assertEquals(HIVDrug.FPV.getDisplayAbbr(), "FPV/r");
		assertEquals(HIVDrug.IDV.getDisplayAbbr(), "IDV/r");
		assertEquals(HIVDrug.LPV.getDisplayAbbr(), "LPV/r");
		assertEquals(HIVDrug.NFV.getDisplayAbbr(), "NFV");
		assertEquals(HIVDrug.SQV.getDisplayAbbr(), "SQV/r");
		assertEquals(HIVDrug.TPV.getDisplayAbbr(), "TPV/r");

		assertEquals(HIVDrug.EFV.getDisplayAbbr(), "EFV");
		assertEquals(HIVDrug.ETR.getDisplayAbbr(), "ETR");
		assertEquals(HIVDrug.NVP.getDisplayAbbr(), "NVP");
		assertEquals(HIVDrug.RPV.getDisplayAbbr(), "RPV");

		assertEquals(HIVDrug.DTG.getDisplayAbbr(), "DTG");
		assertEquals(HIVDrug.EVG.getDisplayAbbr(), "EVG");
		assertEquals(HIVDrug.RAL.getDisplayAbbr(), "RAL");
	}

	@Test
	public void testGetSynonym() {
		assertEquals(HIVDrug.getSynonym("ABC"), HIVDrug.ABC);
		assertEquals(HIVDrug.getSynonym("AZT"), HIVDrug.AZT);
		assertEquals(HIVDrug.getSynonym("D4T"), HIVDrug.D4T);
		assertEquals(HIVDrug.getSynonym("DDI"), HIVDrug.DDI);
		assertEquals(HIVDrug.getSynonym("FTC"), HIVDrug.FTC);
		assertEquals(HIVDrug.getSynonym("3TC"), HIVDrug.LMV);
		assertEquals(HIVDrug.getSynonym("LMV"), HIVDrug.LMV);
		assertEquals(HIVDrug.getSynonym("TDF"), HIVDrug.TDF);

		assertEquals(HIVDrug.getSynonym("ATV/r"), HIVDrug.ATV);
		assertEquals(HIVDrug.getSynonym("DRV/r"), HIVDrug.DRV);
		assertEquals(HIVDrug.getSynonym("FPV/r"), HIVDrug.FPV);
		assertEquals(HIVDrug.getSynonym("IDV/r"), HIVDrug.IDV);
		assertEquals(HIVDrug.getSynonym("LPV/r"), HIVDrug.LPV);
		assertEquals(HIVDrug.getSynonym("ATV"), HIVDrug.ATV);
		assertEquals(HIVDrug.getSynonym("DRV"), HIVDrug.DRV);
		assertEquals(HIVDrug.getSynonym("FPV"), HIVDrug.FPV);
		assertEquals(HIVDrug.getSynonym("IDV"), HIVDrug.IDV);
		assertEquals(HIVDrug.getSynonym("LPV"), HIVDrug.LPV);
		assertEquals(HIVDrug.getSynonym("NFV"), HIVDrug.NFV);
		assertEquals(HIVDrug.getSynonym("SQV/r"), HIVDrug.SQV);
		assertEquals(HIVDrug.getSynonym("TPV/r"), HIVDrug.TPV);
		assertEquals(HIVDrug.getSynonym("SQV"), HIVDrug.SQV);
		assertEquals(HIVDrug.getSynonym("TPV"), HIVDrug.TPV);

		assertEquals(HIVDrug.getSynonym("EFV"), HIVDrug.EFV);
		assertEquals(HIVDrug.getSynonym("ETR"), HIVDrug.ETR);
		assertEquals(HIVDrug.getSynonym("NVP"), HIVDrug.NVP);
		assertEquals(HIVDrug.getSynonym("RPV"), HIVDrug.RPV);

		assertEquals(HIVDrug.getSynonym("DTG"), HIVDrug.DTG);
		assertEquals(HIVDrug.getSynonym("EVG"), HIVDrug.EVG);
		assertEquals(HIVDrug.getSynonym("RAL"), HIVDrug.RAL);
	}

	@Test
	public void testGetSynonymWithException() {
		assertNull(HIVDrug.getSynonym(""));
		assertNull(HIVDrug.getSynonym("EVH"));
	}

	@Test
	public void testValueOf() {
		assertEquals(HIVDrug.valueOf("ABC"), HIVDrug.ABC);
		assertEquals(HIVDrug.valueOf("AZT"), HIVDrug.AZT);
		assertEquals(HIVDrug.valueOf("D4T"), HIVDrug.D4T);
		assertEquals(HIVDrug.valueOf("DDI"), HIVDrug.DDI);
		assertEquals(HIVDrug.valueOf("FTC"), HIVDrug.FTC);
		assertEquals(HIVDrug.valueOf("LMV"), HIVDrug.LMV);
		assertEquals(HIVDrug.valueOf("TDF"), HIVDrug.TDF);

		assertEquals(HIVDrug.valueOf("ATV"), HIVDrug.ATV);
		assertEquals(HIVDrug.valueOf("DRV"), HIVDrug.DRV);
		assertEquals(HIVDrug.valueOf("FPV"), HIVDrug.FPV);
		assertEquals(HIVDrug.valueOf("IDV"), HIVDrug.IDV);
		assertEquals(HIVDrug.valueOf("LPV"), HIVDrug.LPV);
		assertEquals(HIVDrug.valueOf("NFV"), HIVDrug.NFV);
		assertEquals(HIVDrug.valueOf("SQV"), HIVDrug.SQV);
		assertEquals(HIVDrug.valueOf("TPV"), HIVDrug.TPV);

		assertEquals(HIVDrug.valueOf("EFV"), HIVDrug.EFV);
		assertEquals(HIVDrug.valueOf("ETR"), HIVDrug.ETR);
		assertEquals(HIVDrug.valueOf("NVP"), HIVDrug.NVP);
		assertEquals(HIVDrug.valueOf("RPV"), HIVDrug.RPV);

		assertEquals(HIVDrug.valueOf("DTG"), HIVDrug.DTG);
		assertEquals(HIVDrug.valueOf("EVG"), HIVDrug.EVG);
		assertEquals(HIVDrug.valueOf("RAL"), HIVDrug.RAL);
	}
}
