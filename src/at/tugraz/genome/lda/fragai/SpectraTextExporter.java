/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2024 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER. 
 *  
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * by the Free Software Foundation, either version 3 of the License, or 
 * (at your option) any later version.
 *  
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details. 
 *  
 * You should have received a copy of the GNU General Public License 
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 * Please contact lda@genome.tugraz.at if you need additional information or 
 * have any questions.
 */

package at.tugraz.genome.lda.fragai;

import java.io.BufferedOutputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
import java.util.List;

import org.dhatim.fastexcel.Workbook;
import org.dhatim.fastexcel.Worksheet;

import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.vos.SpectrumPointVO;

public class SpectraTextExporter
{
	private final static String HEADER_MZ = "mz Value";
	private final static String HEADER_INTENSITY = "Intensity";
	private ArrayList<SpectrumContainer> spectra_ = null;
	private ArrayList<CombinedSpectrumContainer> combinedSpectra_ = null;
	private String path_;

	public SpectraTextExporter(ArrayList<SpectrumContainer> spectra, String path)
	{
		this.spectra_ = spectra;
		Collections.sort(this.spectra_);
		this.path_ = path;
	}
	
	public SpectraTextExporter(String path, ArrayList<CombinedSpectrumContainer> combinedSpectra)
	{
		this.combinedSpectra_ = combinedSpectra;
		Collections.sort(this.combinedSpectra_);
		this.path_ = path;
	}
	
	public void exportSpectra()
	{
		try (BufferedOutputStream out = new BufferedOutputStream(
				new FileOutputStream(path_));) {
			String s = Settings.VERSION;
			// the constructor can only take a version number in the format xx.yyyy
			Workbook wb = new Workbook(out, "Lipid Data Analyzer",
					s.substring(0, s.indexOf(".", s.indexOf(".") + 1)));
			if (this.spectra_ != null)
			{
				exportRawSpectra(wb);
			}
			else if (this.combinedSpectra_ != null)
			{
				exportCombinedSpectra(wb);
			}
			wb.finish();
			System.out.println("finished");
		}
		catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private void exportCombinedSpectra(Workbook wb)
	{
		for (CombinedSpectrumContainer spectrum : combinedSpectra_)
		{
			ArrayList<SpectrumPointVO> dataPoints = spectrum.computeCombinedSpectrum();
			
			if (!dataPoints.isEmpty())
			{
				String sheetName = String.format("%s_%s", spectrum.getLipidClass(), spectrum.getAdduct());
				Worksheet ws = wb.newWorksheet(sheetName);
				int rowNumber = 0;
				ws.value(rowNumber, 0, "Lipid class:");
				ws.value(rowNumber, 1, spectrum.getLipidClass());
				ws.value(++rowNumber, 0, "Adduct:");
				ws.value(rowNumber, 1, spectrum.getAdduct());
				rowNumber++;
				
				List<String> headerTitles = new ArrayList<String>();
				headerTitles.add(HEADER_MZ);
				headerTitles.add(HEADER_INTENSITY);
				createHeader(ws, headerTitles, ++rowNumber);
				
				
				
				for (SpectrumPointVO dataPoint : dataPoints)
				{
					++rowNumber;
					ws.value(rowNumber, headerTitles.indexOf(HEADER_MZ), dataPoint.getMz());
					ws.value(rowNumber, headerTitles.indexOf(HEADER_INTENSITY), dataPoint.getIntensity());
				}
			}
		}
		System.out.println("done?");
	}
	
	private void exportRawSpectra(Workbook wb)
	{
		int count = 0;
		for (SpectrumContainer spectrum : spectra_)
		{
			Hashtable<Integer,Integer> scanNrLevelHash = spectrum.getScanNrLevelHash();
			for (Integer scanNr : scanNrLevelHash.keySet())
			{
				String sheetName = String.format("%s_%s_%s_%s", spectrum.getEntry().getLipidClass(), spectrum.getAdduct().getAdductName(), scanNrLevelHash.get(scanNr), count++);
				Worksheet ws = wb.newWorksheet(sheetName);
				int rowNumber = 0;
				ws.value(rowNumber, 0, "Lipid class:");
				ws.value(rowNumber, 1, spectrum.getEntry().getLipidClass());
				ws.value(++rowNumber, 0, "Lipid species:");
				ws.value(rowNumber, 1, spectrum.getEntry().getSpecies());
				ws.value(++rowNumber, 0, "Adduct:");
				ws.value(rowNumber, 1, spectrum.getAdduct().getAdductName());
				ws.value(++rowNumber, 0, "Precursor mz:");
				ws.value(rowNumber, 1, spectrum.getScanNrPrecursorHash().get(scanNr).get(0));
				ws.value(++rowNumber, 0, "MS Level:");
				ws.value(rowNumber, 1, scanNrLevelHash.get(scanNr));
				ws.value(++rowNumber, 0, "Scan Number:");
				ws.value(rowNumber, 1, scanNr);
				ws.value(++rowNumber, 0, "Chrom File Name:");
				ws.value(rowNumber, 1, spectrum.getChromFile() == null ? "combined" : spectrum.getChromFile().getName());
				rowNumber++;
				List<String> headerTitles = new ArrayList<String>();
				headerTitles.add(HEADER_MZ);
				headerTitles.add(HEADER_INTENSITY);
				createHeader(ws, headerTitles, ++rowNumber);
				
				ArrayList<SpectrumPointVO> dataPoints = spectrum.getProcessedSpectrum(scanNr);
				
				for (SpectrumPointVO dataPoint : dataPoints)
				{
					++rowNumber;
					ws.value(rowNumber, headerTitles.indexOf(HEADER_MZ), dataPoint.getMz());
					ws.value(rowNumber, headerTitles.indexOf(HEADER_INTENSITY), dataPoint.getIntensity());
				}
			}
		}
	}
	
	
	/**
	 *  Creates a formatted header row with the given header titles in the row
	 * number given by headerRow
	 * 
	 * @param ws
	 *          Excel worksheet to write the header to
	 * @param headerTitles
	 *          List of header titles
	 * @param headerRow
	 * 					the row number for this header
	 */
	private static void createHeader(Worksheet ws, List<String> headerTitles, Integer headerRow)
	{
		for (int i = 0; i < headerTitles.size(); i++) {
			ws.value(headerRow, i, headerTitles.get(i));
		}
		ws.range(headerRow, headerRow, headerRow, headerTitles.size()).style()
				.bold().horizontalAlignment("center").fontName("Arial").fontSize(12)
				.set();
	}
	
	
}
