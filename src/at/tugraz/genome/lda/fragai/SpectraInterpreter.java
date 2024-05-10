package at.tugraz.genome.lda.fragai;

import java.io.File;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Vector;

import at.tugraz.genome.lda.vos.SpectrumPointVO;
import at.tugraz.genome.maspectras.quantification.CgProbe;

public class SpectraInterpreter
{
	private ArrayList<SpectrumContainer> spectra_;
	
	public SpectraInterpreter(ArrayList<SpectrumContainer> spectra)
	{
		this.spectra_ = spectra;
	}
	
	/**
	 * TODO: only take the highest intensity MZ for each accumulation of peaks, (and add all other intensities of this accumulation to it!?)
	 * After this step we only have the peakmax plus total intensity.
	 * Get a variable that gets us the deviation for each of these accumulations for later merging.
	 * This is then used to compare spectra!
	 * We generate cumulative spectra with same Class and Adduct and output the fragments (for now) found in all cases.
	 * Then add the precursorcleared option (removing the precursor for all fragments) and creating cumulative spectra this way (adding this to the other constant fragments).
	 * 
	 * 
	 */
	public ArrayList<SpectrumContainer> interpretSpectra()
	{
		Hashtable<String, ArrayList<SpectrumContainer>> spectraForLipid = groupSpectraContainer();
		ArrayList<CombinedSpectrumContainer> combinedSpectra = new ArrayList<CombinedSpectrumContainer>();
		for (String lipid : spectraForLipid.keySet())
		{
			CombinedSpectrumContainer combined = new CombinedSpectrumContainer(lipid.substring(0, lipid.indexOf("_")), lipid.substring(lipid.indexOf("_")+1, lipid.length()));
			for (SpectrumContainer container : spectraForLipid.get(lipid))
			{
				combined.addContainerCopy(container);
			}
			
			System.out.println("HI!");
		}
		
		
		System.out.println("HI");
		return spectra_;
	}
	
	private String getSpectrumID(SpectrumContainer spectrum)
	{
		return spectrum.getEntry().getLipidClass()+"_"+spectrum.getAdduct().getAdductName();
	}
	
	private Hashtable<String, ArrayList<SpectrumContainer>> groupSpectraContainer()
	{
		Hashtable<String, ArrayList<SpectrumContainer>> spectraForLipid = new Hashtable<String, ArrayList<SpectrumContainer>>();
		for (SpectrumContainer spectrum : spectra_)
		{
			String id = getSpectrumID(spectrum);
			if (!spectraForLipid.containsKey(id))
			{
				spectraForLipid.put(id, new ArrayList<SpectrumContainer>());
			}
			spectraForLipid.get(id).add(spectrum);
		}
		return spectraForLipid;
	}
	
	/**
	 * go through spectrum containers. group them according to lipid class and adduct.
	 * 
	 * this should contain: entry, adduct
	 * list of spectrum containers
	 * 
	 * 
	 * for each theoretical 
	 * a list of SpectrumPointVOs 
	 * 
	 * 
	 */
	private class CombinedSpectrumContainer
	{
		private String lipidClass_;
		private String adduct_;
		private ArrayList<SpectrumContainer> containers_;
		
		private CombinedSpectrumContainer(String lipidClass, String adduct)
		{
			this.lipidClass_ = lipidClass;
			this.adduct_ = adduct;
			this.containers_ = new ArrayList<SpectrumContainer>();
		}
		
		/**
		 * Checks if there are scans of at least 2 different precursor masses are available for ms level == 2.
		 * There must be a total of at least 3 scans.
		 * @return true if the criteria are fulfilled.
		 */
		private boolean isViableData()
		{
//			Hashtable<String,Integer> formulas = 
			for (SpectrumContainer container : containers_)
			{
				container.getEntry().getSumFormula();
			}
			return true;
		}
		
		
		
		
		/**
		 * Adds a copy of the container clearing the spectra data up. Only the peak maxima are saved.
		 * For now: only MS2 scans are used.
		 * @param container
		 */
		private void addContainerCopy(SpectrumContainer container)
		{
			if (container.getEntry().getLipidClass().equalsIgnoreCase(lipidClass_) && container.getAdduct().getAdductName().equalsIgnoreCase(adduct_))
			{
				SpectrumContainer newContainer = new SpectrumContainer(container);
				for (Integer scanNr : newContainer.getScanNrLevelHash().keySet())
				{
					ArrayList<Integer> scans = new ArrayList<Integer>();
					scans.add(scanNr);
					ArrayList<SpectrumPointVO> dataPoints = clearData(newContainer.getProcessedSpectrum(scanNr));
					newContainer.setProcessedSpectrum(scanNr, dataPoints);
				}
				this.containers_.add(newContainer);
			}
		}
		
		private ArrayList<SpectrumPointVO> clearData(ArrayList<SpectrumPointVO> dataPoints)
		{
			ArrayList<SpectrumPointVO> clearedData = new ArrayList<SpectrumPointVO>();
			SpectrumPointVO max = new SpectrumPointVO("", 0.0f);
			for (SpectrumPointVO vo : dataPoints)
			{
				if (vo.getIntensity()>0)
				{
					if (vo.getIntensity() > max.getIntensity())
					{
						max = vo;
					}
				}
				else if (max.getIntensity()>0)
				{
					clearedData.add(max);
					max = new SpectrumPointVO("", 0.0f);
				}
			}
			return clearedData;
		}
		
	}
	
}
