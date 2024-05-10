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
	private ArrayList<SpectrumContainer> spectra_;
	private String path_;

	public SpectraTextExporter(ArrayList<SpectrumContainer> spectra, String path)
	{
		this.spectra_ = spectra;
		Collections.sort(this.spectra_);
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
					ws.value(rowNumber, 1, spectrum.getChromFile().getName());
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
