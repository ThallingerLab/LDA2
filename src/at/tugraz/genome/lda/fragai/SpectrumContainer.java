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

import java.io.File;
import java.nio.ByteBuffer;
import java.nio.FloatBuffer;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Hashtable;
import java.util.Objects;
import java.util.Vector;

import at.tugraz.genome.dbutilities.Base64;
import at.tugraz.genome.lda.vos.AdductVO;
import at.tugraz.genome.lda.vos.SpectrumPointVO;
import at.tugraz.genome.maspectras.quantification.CgProbe;

public class SpectrumContainer implements Comparable<SpectrumContainer>
{
	FragTargetListEntry entry_;
	AdductVO adduct_;
	CgProbe probe_;
	File chromFile_;
	/** key is the scan number, value the spectrum in Base64 encoding */
	Hashtable<Integer,String> scanNrSpectrumHash_;
	/** key is the scan number, value the processed spectrum */
	Hashtable<Integer,ArrayList<SpectrumPointVO>> scanNrProcessedSpectrumHash_;
	/** key is the scan number, value the precursor m/z value */
	Hashtable<Integer,Vector<Double>> scanNrPrecursorHash_;
	/** key is the scan number, value the ms level */
	Hashtable<Integer,Integer> scanNrLevelHash_;
	
	public SpectrumContainer(FragTargetListEntry entry, AdductVO adduct,
			CgProbe probe, File chromFile, Hashtable<Integer,String> scanNrSpectrumHash,
			Hashtable<Integer,Vector<Double>> scanNrPrecursorHash,
			Hashtable<Integer,Integer> scanNrLevelHash)
	{
		this.entry_ = entry;
		this.adduct_ = adduct;
		this.probe_ = probe;
		this.chromFile_ = chromFile;
		this.scanNrSpectrumHash_ = scanNrSpectrumHash;
		this.scanNrPrecursorHash_ = scanNrPrecursorHash;
		this.scanNrLevelHash_ = scanNrLevelHash;
		this.scanNrProcessedSpectrumHash_ = processSpectra();
	}
	
	public SpectrumContainer(SpectrumContainer other)
	{
		this(other.getEntry(), other.getAdduct(), other.getProbe(), other.getChromFile(), other.getScanNrSpectrumHash(), other.getScanNrPrecursorHash(), other.getScanNrLevelHash());
	}
	
	private Hashtable<Integer,ArrayList<SpectrumPointVO>> processSpectra()
	{
		Hashtable<Integer,ArrayList<SpectrumPointVO>> processedSpectra = new Hashtable<Integer,ArrayList<SpectrumPointVO>>();
		for (Integer scanNr : getScanNrLevelHash().keySet())
		{
			ArrayList<Integer> scans = new ArrayList<Integer>();
			scans.add(scanNr);
			ArrayList<SpectrumPointVO> dataPoints = mergeSpectra(scans);
			processedSpectra.put(scanNr, dataPoints);
		}
		return processedSpectra;
	}
	
	private ArrayList<SpectrumPointVO> mergeSpectra(Collection<Integer> scanNumbers)
	{
    Hashtable<String,SpectrumPointVO> intensityPoints = new Hashtable<String,SpectrumPointVO>();  
    for (Integer  scanNr : scanNumbers){
      String spectrum = scanNrSpectrumHash_.get(scanNr);
      FloatBuffer buffer = ByteBuffer.wrap(Base64.decode(spectrum)).asFloatBuffer();
      int iItems = buffer.limit();
      float mzValue = 0;
      for(int iItem = 0; iItem < iItems; iItem++){
        if (iItem%2==0) mzValue = buffer.get();
        else{
          float intensity = buffer.get();
          String mzOriginal = String.valueOf(mzValue);
          SpectrumPointVO vo = null;
          if (intensityPoints.containsKey(mzOriginal))
            vo = intensityPoints.get(mzOriginal);
          else
            vo = new SpectrumPointVO(mzOriginal,mzValue);
          vo.addIntensity(intensity);
          intensityPoints.put(mzOriginal, vo);
        }
      }
    }  
    ArrayList<SpectrumPointVO> vos = new ArrayList<SpectrumPointVO>(intensityPoints.values());
    Collections.sort(vos);
    return vos;
  }
	
	

	public FragTargetListEntry getEntry()
	{
		return entry_;
	}


	public AdductVO getAdduct()
	{
		return adduct_;
	}


	public CgProbe getProbe()
	{
		return probe_;
	}


	public File getChromFile()
	{
		return chromFile_;
	}


	public Hashtable<Integer,String> getScanNrSpectrumHash()
	{
		return scanNrSpectrumHash_;
	}

	public Hashtable<Integer,Vector<Double>> getScanNrPrecursorHash()
	{
		return scanNrPrecursorHash_;
	}


	public Hashtable<Integer,Integer> getScanNrLevelHash()
	{
		return scanNrLevelHash_;
	}
	

	public ArrayList<SpectrumPointVO> getProcessedSpectrum(Integer scanNumber)
	{
		return scanNrProcessedSpectrumHash_.get(scanNumber);
	}
	
	public void setProcessedSpectrum(Integer scanNumber, ArrayList<SpectrumPointVO> spectrum)
	{
		scanNrProcessedSpectrumHash_.put(scanNumber, spectrum);
	}
	
	

	@Override
	public int hashCode()
	{
		return Objects.hash(adduct_, chromFile_, entry_, probe_, scanNrLevelHash_,
				scanNrPrecursorHash_, scanNrSpectrumHash_);
	}

	@Override
	public boolean equals(Object obj)
	{
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		SpectrumContainer other = (SpectrumContainer) obj;
		return Objects.equals(adduct_, other.adduct_)
				&& Objects.equals(chromFile_, other.chromFile_)
				&& Objects.equals(entry_, other.entry_)
				&& Objects.equals(probe_, other.probe_)
				&& Objects.equals(scanNrLevelHash_, other.scanNrLevelHash_)
				&& Objects.equals(scanNrPrecursorHash_, other.scanNrPrecursorHash_)
				&& Objects.equals(scanNrSpectrumHash_, other.scanNrSpectrumHash_);
	}
	
  public int compareTo(SpectrumContainer other) 
  {
  	return Comparator
  			.comparing((SpectrumContainer sc) -> sc.getEntry().getLipidClass())
  			.thenComparing((SpectrumContainer sc) -> sc.getEntry().getSpecies())
  			.thenComparing((SpectrumContainer sc) -> sc.getEntry().hashCode())
  			.compare(this, other);
  }
	
}
