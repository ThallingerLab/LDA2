package at.tugraz.genome.lda.fragai;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;

import at.tugraz.genome.lda.vos.SpectrumPointVO;

/**
 * Additionally to the default data of a spectrum point, this class contains information of 
 * how often a point was detected in a group of spectra
 * the average intensity
 * the set of lipid species to determine if this spectrum point is supported by statistics
 * 
 * @author Leonida M. Lamp
 */
public class CombinedSpectrumPointVO implements Comparable<CombinedSpectrumPointVO>
{
	private Boolean precursorAdjusted_;
	private HashMap<String,Integer> speciesSpectraCountMap_;
	private HashMap<String,ArrayList<Double>> speciesIntensityMap_;
	private ArrayList<Double> precursorMZ_;

	public CombinedSpectrumPointVO(Boolean precursorRemoved)
	{
		precursorAdjusted_ = precursorRemoved;
		speciesSpectraCountMap_ = new HashMap<String,Integer>();
		speciesIntensityMap_ = new HashMap<String,ArrayList<Double>>();
		precursorMZ_ = new ArrayList<Double>();
	}
	
	public void addSpectrumPoint(SpectrumPointVO vo, String species)
	{
		addToSpeciesSpectraCountMap(species);
		addToSpeciesIntensityMap(species, vo.getIntensity());
		addToPrecursorMZ(vo.getMz());
	}
	
	private void addToSpeciesSpectraCountMap(String species)
	{
		if (!speciesSpectraCountMap_.containsKey(species))
		{
			speciesSpectraCountMap_.put(species, 1);
		}
		else
		{
			int num = speciesSpectraCountMap_.get(species);
			speciesSpectraCountMap_.put(species, ++num);
		}
	}
	
	private void addToSpeciesIntensityMap(String species, float intensity)
	{
		if (!speciesIntensityMap_.containsKey(species))
		{
			speciesIntensityMap_.put(species, new ArrayList<Double>());
		}
		speciesIntensityMap_.get(species).add(new Double(intensity));
	}
	
	private void addToPrecursorMZ(float value)
	{
		precursorMZ_.add(new Double(value));
	}
	
	public double getAverageIntensity()
	{
		double sum = 0.0;
		int count = 0;
		for (String species : speciesIntensityMap_.keySet())
		{
			for (Double intensity : speciesIntensityMap_.get(species))
			{
				sum += intensity; 
				count++;
			}
		}
		return sum/count;
	}
	
	public double getAverageMZ()
	{
		double sum = 0.0;
		for (Double mz : precursorMZ_)
		{
			sum += mz;
		}
		return sum/precursorMZ_.size();
	}
	
	public Integer getTotalCount()
	{
		int sum = 0;
		for (String species : speciesSpectraCountMap_.keySet())
		{
			sum += speciesSpectraCountMap_.get(species);
		}
		return sum;
	}
	
	public int getSpeciesNumber()
	{
		return speciesSpectraCountMap_.keySet().size();
	}
	
	public String getPerSpeciesCount()
	{
		StringBuilder builder = new StringBuilder();
		for (String species : speciesSpectraCountMap_.keySet())
		{
			builder.append(species+": "+speciesSpectraCountMap_.get(species)+"; ");
		}
		return builder.toString();
	}
	
	public boolean getPrecursorAdjusted()
	{
		return this.precursorAdjusted_;
	}

	@Override
	public int compareTo(CombinedSpectrumPointVO anotherVO)
	{
		return Comparator
				.comparing((CombinedSpectrumPointVO vo) -> vo.getSpeciesNumber()).reversed()
				.thenComparing((CombinedSpectrumPointVO vo) -> vo.getPrecursorAdjusted())
  			.thenComparing((CombinedSpectrumPointVO vo) -> vo.getAverageMZ())
  			.compare(this,anotherVO);
	}

}
