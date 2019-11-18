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

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.*;

import org.junit.Test;

public class DrugClassTest {

	@Test
	public void testGetGene() {
		assertEquals(HIVDrugClass.PI.gene(), HIVGeneEnum.PR);
		assertEquals(HIVDrugClass.NRTI.gene(), HIVGeneEnum.RT);
		assertEquals(HIVDrugClass.NNRTI.gene(), HIVGeneEnum.RT);
		assertEquals(HIVDrugClass.INSTI.gene(), HIVGeneEnum.IN);
	}

	@Test
	public void testGetSynonym() {
		assertEquals(HIVDrugClass.getSynonym("PI"), HIVDrugClass.PI);
		assertEquals(HIVDrugClass.getSynonym("NRTI"), HIVDrugClass.NRTI);
		assertEquals(HIVDrugClass.getSynonym("NNRTI"), HIVDrugClass.NNRTI);
		assertEquals(HIVDrugClass.getSynonym("INSTI"), HIVDrugClass.INSTI);
		assertEquals(HIVDrugClass.getSynonym("INI"), HIVDrugClass.INSTI);
	}

	@Test
	public void testGetAllDrugs() {
		List<HIVDrug> piExpecteds = new ArrayList<>();
		piExpecteds.add(HIVDrug.ATV);
		piExpecteds.add(HIVDrug.DRV);
		piExpecteds.add(HIVDrug.FPV);
		piExpecteds.add(HIVDrug.IDV);
		piExpecteds.add(HIVDrug.LPV);
		piExpecteds.add(HIVDrug.NFV);
		piExpecteds.add(HIVDrug.SQV);
		piExpecteds.add(HIVDrug.TPV);
		assertEquals(
			HIVDrugClass.PI.getAllDrugs(), piExpecteds);

		List<HIVDrug> nrtiExpecteds = new ArrayList<>();
		nrtiExpecteds.add(HIVDrug.ABC);
		nrtiExpecteds.add(HIVDrug.AZT);
		nrtiExpecteds.add(HIVDrug.D4T);
		nrtiExpecteds.add(HIVDrug.DDI);
		nrtiExpecteds.add(HIVDrug.FTC);
		nrtiExpecteds.add(HIVDrug.LMV);
		nrtiExpecteds.add(HIVDrug.TDF);
		assertEquals(
			HIVDrugClass.NRTI.getAllDrugs(), nrtiExpecteds);

		List<HIVDrug> nnrtiExpecteds = new ArrayList<>();
		nnrtiExpecteds.add(HIVDrug.DOR);
		nnrtiExpecteds.add(HIVDrug.EFV);
		nnrtiExpecteds.add(HIVDrug.ETR);
		nnrtiExpecteds.add(HIVDrug.NVP);
		nnrtiExpecteds.add(HIVDrug.RPV);
		assertEquals(
			HIVDrugClass.NNRTI.getAllDrugs(), nnrtiExpecteds);

		List<HIVDrug> instiExpecteds = new ArrayList<>();
		instiExpecteds.add(HIVDrug.BIC);
		instiExpecteds.add(HIVDrug.DTG);
		instiExpecteds.add(HIVDrug.EVG);
		instiExpecteds.add(HIVDrug.RAL);
		assertEquals(
			HIVDrugClass.INSTI.getAllDrugs(), instiExpecteds);
	}

	@Test
	public void testGetFullName() {
		assertEquals(
			"Nucleoside Reverse Transcriptase Inhibitor",
			HIVDrugClass.NRTI.getFullName());
		assertEquals(
			"Non-nucleoside Reverse Transcriptase Inhibitor",
			HIVDrugClass.NNRTI.getFullName());
		assertEquals(
			"Protease Inhibitor", HIVDrugClass.PI.getFullName());
		assertEquals(
			"Integrase Strand Transfer Inhibitor",
			HIVDrugClass.INSTI.getFullName());
	}
}
