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

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;

import org.apache.commons.math3.util.Pair;

import at.tugraz.genome.lda.vos.SpectrumPointVO;


/**
 * 
 * @author Leonida M. Lamp
 *
 */
public class CombinedSpectrumContainer implements Comparable<CombinedSpectrumContainer>
{
	private String lipidClass_;
	private String adduct_;
	private Integer mzTolerance_;
	private ArrayList<SpectrumContainer> containers_;
	
	public CombinedSpectrumContainer(String lipidClass, String adduct, Integer mzTolerance)
	{
		this.lipidClass_ = lipidClass;
		this.adduct_ = adduct;
		this.mzTolerance_ = mzTolerance;
		this.containers_ = new ArrayList<SpectrumContainer>();
	}
	
	/**
	 * Checks if there are scans of at least 2 different precursor masses are available for ms level == 2.
	 * There must be a total of at least 3 scans.
	 * @return true if the criteria are fulfilled.
	 */
	protected boolean isViableData()
	{
		HashSet<Double> uniqueMolMasses = new HashSet<Double>();
		Integer count = 0;
		for (SpectrumContainer container : containers_)
		{
			uniqueMolMasses.add(container.getEntry().computeTheoreticalPrecursorMZValue(adduct_));
			for (Integer level : container.getScanNrLevelHash().values())
			{
				if (level == 2) count ++;
			}
		}
		if (count>2 && uniqueMolMasses.size()>1)
			return true;
		return false;
	}
	
	/**
	 * Adds a copy of the container clearing the spectra data up. Only the peak maxima for each accumulation of intensities are saved.
	 * For now: only MS2 scans are used.
	 * @param container
	 */
	protected void addContainerCopy(SpectrumContainer container)
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
	
	/**
	 * Returns a cumulative spectrum for this container.
	 * First, all data points in common are added to the cumulative spectrum 
	 * With intensities relative to the overall maximum intensity for each scan. (Anther idea reasonable?)
	 * Then, the precursorcleared option is added (the precursor is removed for all fragements)
	 * @return sorted datapoints of the cumulative spectrum.
	 */
	protected ArrayList<SpectrumPointVO> computeCombinedSpectrum()
	{
		ArrayList<SpectrumPointVO> combined = new ArrayList<SpectrumPointVO>();
		HashMap<Pair<Double,Double>,ArrayList<SpectrumPointVO>> groups = new HashMap<Pair<Double,Double>,ArrayList<SpectrumPointVO>>();
		int countScan = 0;
		for (SpectrumContainer container : containers_)
		{
			for (Integer scanNr : container.getScanNrLevelHash().keySet())
			{
				if (container.getScanNrLevelHash().get(scanNr) == 2)
				{
					countScan++;
					ArrayList<SpectrumPointVO> processedSpectrum = container.getProcessedSpectrum(scanNr);
					float highestInt = 0f;
					for (SpectrumPointVO point : processedSpectrum)
					{
						highestInt = highestInt < point.getIntensity() ? point.getIntensity() : highestInt;
					}
					
					boolean added = false;
					
					for (SpectrumPointVO point : processedSpectrum)
					{
						ArrayList<Pair<Double, Double>> bins = new ArrayList<Pair<Double, Double>>(groups.keySet());
						for (int i=0; i<bins.size(); i++)
						{
							Pair<Double,Double> bin = bins.get(i);
							if (bin.getFirst() <= point.getMz() && bin.getSecond() >= point.getMz())
							{
								addPointToBin(bin, groups, point, highestInt);
								added = true;
							}
						}
						if (!added)
						{
							addPointToBin(null, groups, point, highestInt);
						}
					}
				}
			}
		}
		System.out.println(countScan);System.out.println(countScan/2+1);
		for (Pair<Double,Double> bin : groups.keySet())
		{
			ArrayList<SpectrumPointVO> points = groups.get(bin);
			if (points.size() > countScan/2+1)
			{
				float averageMZ = 0f;
				float averageIntensity = 0f;
				
				for (SpectrumPointVO point : points)
				{
					averageMZ += point.getMz();
					averageIntensity += point.getIntensity();
				}
				averageMZ /= points.size();
				averageIntensity /= points.size();
				
				SpectrumPointVO averaged = new SpectrumPointVO(averageMZ, averageIntensity);
				combined.add(averaged);
			}
		}
		Collections.sort(combined);
		
		ArrayList<SpectrumPointVO> combinedTemp = new ArrayList<SpectrumPointVO>();
		HashMap<Pair<Double,Double>,ArrayList<SpectrumPointVO>> groupsTemp = new HashMap<Pair<Double,Double>,ArrayList<SpectrumPointVO>>();
		for (SpectrumContainer container : containers_)
		{
			float precursorMz = new Float(container.getEntry().computeTheoreticalPrecursorMZValue(adduct_));
			for (Integer scanNr : container.getScanNrLevelHash().keySet())
			{
				if (container.getScanNrLevelHash().get(scanNr) == 2)
				{
					ArrayList<SpectrumPointVO> processedSpectrum = container.getProcessedSpectrum(scanNr);
					float highestInt = 0f;
					for (SpectrumPointVO point : processedSpectrum)
					{
						highestInt = highestInt < point.getIntensity() ? point.getIntensity() : highestInt;
					}
					
					boolean added = false;
					
					for (SpectrumPointVO pointUnadjusted : processedSpectrum)
					{
						SpectrumPointVO point = new SpectrumPointVO(precursorMz-pointUnadjusted.getMz(), pointUnadjusted.getIntensity());
						ArrayList<Pair<Double, Double>> bins = new ArrayList<Pair<Double, Double>>(groupsTemp.keySet());
						for (int i=0; i<bins.size(); i++)
						{
							Pair<Double,Double> bin = bins.get(i);
							if (bin.getFirst() <= point.getMz() && bin.getSecond() >= point.getMz())
							{
								addPointToBin(bin, groupsTemp, point, highestInt);
								added = true;
							}
						}
						if (!added)
						{
							addPointToBin(null, groupsTemp, point, highestInt);
						}
					}
				}
			}
		}
		System.out.println(countScan);System.out.println(countScan/2+1);
		for (Pair<Double,Double> bin : groupsTemp.keySet())
		{
			ArrayList<SpectrumPointVO> points = groupsTemp.get(bin);
			if (points.size() > countScan/2+1)
			{
				float averageMZ = 0f;
				float averageIntensity = 0f;
				
				for (SpectrumPointVO point : points)
				{
					averageMZ += point.getMz();
					averageIntensity += point.getIntensity();
				}
				averageMZ /= points.size();
				averageIntensity /= points.size();
				
				SpectrumPointVO averaged = new SpectrumPointVO(averageMZ, averageIntensity);
				combinedTemp.add(averaged);
			}
		}
		Collections.sort(combinedTemp);
		combined.addAll(combinedTemp);
		
		return combined;
	}
	
	private void addPointToBin(Pair<Double,Double> bin, HashMap<Pair<Double,Double>,ArrayList<SpectrumPointVO>> groups, SpectrumPointVO newPoint, float highestInt)
	{
		SpectrumPointVO adjustedInt = new SpectrumPointVO(newPoint.getMz(), newPoint.getIntensity()/highestInt);
		if (bin != null)
		{
			ArrayList<SpectrumPointVO> points = groups.get(bin);
			groups.remove(bin);
			points.add(adjustedInt);
			float sumMz = 0f;
			for (SpectrumPointVO point : points)
			{
				sumMz += point.getMz();
			}
			float averageMz = sumMz/points.size();
			Pair<Double,Double> newBin = computeLowerUpperPair(averageMz);
			groups.put(newBin, points);
		}
		else
		{
			ArrayList<SpectrumPointVO> points = new ArrayList<SpectrumPointVO>();
			points.add(adjustedInt);
			Pair<Double,Double> newBin = computeLowerUpperPair(adjustedInt.getMz());
			groups.put(newBin, points);
		}
	}
	
	private Pair<Double,Double> computeLowerUpperPair(float mz)
	{
		Double tolerance = computeTolerance(mz);
		Double lowerLimit = mz-tolerance;
		Double upperLimit = mz+tolerance;
		return new Pair<Double,Double>(lowerLimit,upperLimit);
	}
	
	private Double computeTolerance(float mz)
	{
		return (double) (mz/1000000*mzTolerance_);
	}

	public String getLipidClass()
	{
		return lipidClass_;
	}

	public String getAdduct()
	{
		return adduct_;
	}
	

	@Override
	public int compareTo(CombinedSpectrumContainer other)
	{
		return Comparator
  			.comparing(CombinedSpectrumContainer::getLipidClass)
  			.thenComparing(CombinedSpectrumContainer::getAdduct)
  			.compare(this, other);
	}
	
}
