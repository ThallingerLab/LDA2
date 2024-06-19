/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2023 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp
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

package at.tugraz.genome.lda.target.experiment;

import java.util.HashSet;
import java.util.Hashtable;
import java.util.Set;
import java.util.Vector;

import javax.swing.JFrame;

import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.quantification.QuantificationResult;
import at.tugraz.genome.lda.target.JDefaultComponents;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.ResultFileVO;
import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
import at.tugraz.genome.lda.target.IsotopeLabelVO;

/**
 * 
 * @author Leonida M. Lamp
 *
 */
public class ExperimentLabelDefinitionPanel extends ExperimentTableInputPanel
{
	private static final long serialVersionUID = 1L;
	
	private static final String ELEMENT_D = "D";
	private static final String ELEMENT_13C = "Cc";
	private static final int EDITABLE_COLUMN = 3;
	
	private Vector<IsotopeLabelVO> unambiguousIsotopeLabels_;
  
  public ExperimentLabelDefinitionPanel(JDefaultComponents wizardComponents, ExperimentResultFileChooserPanel experimentFileChooserPanel) {
      super(wizardComponents, experimentFileChooserPanel, 
      		"Create RT-DB using experimental data of SIL specific for \u03C9-positions.", "Enter the \u03C9-position associated with each detected SIL.");
  }
  
  @Override
  public void initDataDisplay()
  {
  	cleanPanels();
  	String[] columnNames = { "Label ID", "Isotope", "Amount of isotope atoms", "Enter \u03C9-position here!"};
    Object[][] tableData = generateTableData();
  	getDefaultComponents().updateComponents();
  	init(generateDisplayPanel(columnNames, tableData, EDITABLE_COLUMN));
  }
  
  @Override
  protected Object[][] generateTableData()
  {
  	Object[][] tableData = new Object[unambiguousIsotopeLabels_.size()][4];
  	
  	int count=0;
    for (IsotopeLabelVO labelVO : unambiguousIsotopeLabels_)
    {
    	Hashtable<String,Integer> labelElements = labelVO.getLabelElements();
    	String isotope;
    	int isotopeAmount;
    	if (labelElements.containsKey(ELEMENT_D) && labelElements.containsKey(ELEMENT_13C))
    	{
    		isotope = "<html> D and <sup>13</sup>C </html>";
    		isotopeAmount = labelElements.get(ELEMENT_D) + labelElements.get(ELEMENT_13C);
    	}
    	else if (labelElements.containsKey(ELEMENT_D))
    	{
    		isotope = "D";
    		isotopeAmount = labelElements.get(ELEMENT_D);
    	}
    	else
    	{
    		isotope = "<html> <sup>13</sup>C </html>";
    		isotopeAmount = labelElements.get(ELEMENT_13C);
    	}
    	tableData[count][0] = labelVO.getLabelId();
    	tableData[count][1] = isotope;
    	tableData[count][2] = isotopeAmount;
    	tableData[count][EDITABLE_COLUMN] = 0;
    	count++;
    }
    
    return tableData;
  }
  
  protected void parseDataForLabels()
  {
  	Hashtable<String, IsotopeLabelVO> isotopeLabels = new Hashtable<String, IsotopeLabelVO>();
  	Set<String> ambiguousLabels = new HashSet<String>();
  	for (ResultFileVO resultFileVO : getResultFiles())
  	{
  		QuantificationResult quantRes = resultFileVO.getQuantificationResult();
  		Hashtable<String, Vector<LipidParameterSet>> identifications = quantRes.getIdentifications();
  		for (String key : identifications.keySet())
			{
				for (LipidParameterSet analyte : identifications.get(key))
				{
					String analyteName = analyte.getName();
					if (containsPotentialSingleLetterLabel(analyteName))
					{
						String labelID = analyteName.substring(0, 1);
						IsotopeLabelVO labelVO;
						Hashtable<String,Integer> categorized = null;
						try {categorized = StaticUtils.categorizeFormula(analyte.getAnalyteFormula());} catch (Exception ex) {}
						
						if (isotopeLabels.containsKey(labelID))
						{
							labelVO = isotopeLabels.get(labelID);
							if (!isLabelIdentical(labelVO, categorized))
							{
								ambiguousLabels.add(labelID);
							}
						}
						else
						{
							Hashtable<String,Integer> labelElements = initLabelElements(categorized);
							labelVO = new IsotopeLabelVO(labelElements, -1, labelID);
							isotopeLabels.put(labelID, labelVO);
						}
					}
				}
			}
  	}
  	unambiguousIsotopeLabels_ = filterAmbiguousLabels(isotopeLabels, ambiguousLabels);
  	if (unambiguousIsotopeLabels_.size() < 1)
  	{
  		new WarningMessage(new JFrame(), "Warning", "<html> No unambiguous stable isotope labels were detected. <br>"
  				+ "Label IDs must be uppercase letters (A to Z) and correspond to a defined number of D and/or <sup>13</sup>C elements. </html>");
			goBack();
  	}
  }
  
  private Vector<IsotopeLabelVO> filterAmbiguousLabels(Hashtable<String, IsotopeLabelVO> isotopeLabels, Set<String> ambiguousLabels)
  {
  	Vector<IsotopeLabelVO> unambiguousLabels = new Vector<IsotopeLabelVO>();
  	for (String labelID : isotopeLabels.keySet())
  	{
  		if (!ambiguousLabels.contains(labelID))
  		{
  			unambiguousLabels.add(isotopeLabels.get(labelID));
  		}
  	}
  	return unambiguousLabels;
  }
  
  private Hashtable<String,Integer> initLabelElements(Hashtable<String,Integer> categorized)
  {
  	Hashtable<String,Integer> labelElements = new Hashtable<String,Integer>();
		if (categorized.keySet().contains(ELEMENT_D))
		{
			labelElements.put("H", -categorized.get(ELEMENT_D));
			labelElements.put(ELEMENT_D, categorized.get(ELEMENT_D));
		}
		if (categorized.keySet().contains(ELEMENT_13C))
		{
			labelElements.put("C", -categorized.get(ELEMENT_13C));
			labelElements.put(ELEMENT_13C, categorized.get(ELEMENT_13C));
		}
		return labelElements;
  }
  
  private boolean isLabelIdentical(IsotopeLabelVO labelVO, Hashtable<String,Integer> categorized)
  {
  	Hashtable<String,Integer> labelElements = labelVO.getLabelElements();
  	if (labelElements.keySet().contains(ELEMENT_D) && categorized.keySet().contains(ELEMENT_D))
  	{
  		if (labelElements.get(ELEMENT_D) == categorized.get(ELEMENT_D))
  		{
  			return true;
  		}
  	}
  	if (labelElements.keySet().contains(ELEMENT_13C) && categorized.keySet().contains(ELEMENT_13C))
  	{
  		if (labelElements.get(ELEMENT_13C) == categorized.get(ELEMENT_13C))
  		{
  			return true;
  		}
  	}
  	return false;
  }
  
  private boolean containsPotentialSingleLetterLabel(String analyteName)
  {
  	if (analyteName.length()>1 && isUpperCaseLetter(analyteName.charAt(0)) && !isUpperCaseLetter(analyteName.charAt(1)))
  	{
  		return true;
  	}
  	return false;
  }
  
  private boolean isUpperCaseLetter(char c)
  {
  	if (c >= 'A' && c <= 'Z') 
  	{
  		return true;
  	}
  	return false;
  }
  
  //TODO: this is not part of the main code and should be deleted when not needed anymore
  protected void parseDataForRTCalibrationAnchors()
  {
  	Hashtable<String, Vector<CalibrationAnchor>> molecularSpeciesPerResult = new Hashtable<String, Vector<CalibrationAnchor>>();
  	for (ResultFileVO resultFileVO : getResultFiles())
  	{
  		QuantificationResult quantRes = resultFileVO.getQuantificationResult();
  		Hashtable<String, Vector<LipidParameterSet>> identifications = quantRes.getIdentifications();
  		for (String key : identifications.keySet())
			{
				for (LipidParameterSet analyte : identifications.get(key))
				{
					if (analyte instanceof LipidomicsMSnSet)
					{
						LipidomicsMSnSet analyteMSn = (LipidomicsMSnSet) analyte;
						Set<String> molecularSpeciesSet = analyteMSn.getPositionInsensitiveHumanReadableNameSet(false);
						if (molecularSpeciesSet.size() == 1) //we only want peaks containing just one molecular species
						{
							String molecularSpecies = key+" "+molecularSpeciesSet.iterator().next();	
							CalibrationAnchor anchor = new CalibrationAnchor(molecularSpecies, resultFileVO.getFileName(), analyteMSn);
							if (!molecularSpeciesPerResult.containsKey(molecularSpecies))
							{
								molecularSpeciesPerResult.put(molecularSpecies, new Vector<CalibrationAnchor>());
							}
							molecularSpeciesPerResult.get(molecularSpecies).add(anchor);
						}
					}
				}
			}
  	}
  	System.out.println("!!! new !!!");
  	int count = 0;
  	for (String molecularSpecies : molecularSpeciesPerResult.keySet())
  	{
  		Vector<CalibrationAnchor> anchors = molecularSpeciesPerResult.get(molecularSpecies);
  		if (anchors.size() == getResultFiles().size()) //we only want molecular species found in all files
  		{
  			Set<String> uniqueFileNames = new HashSet<String>();
  			for (CalibrationAnchor anchor : anchors)
  			{
  				uniqueFileNames.add(anchor.getFilePath());
  			}
  			if (uniqueFileNames.size() == getResultFiles().size()) //we only want molecular species only found once per file
  			{
  				for (CalibrationAnchor anchor : anchors)
    			{
    				System.out.println(molecularSpecies+" \t File: "+anchor.getFilePath()+" \t RT: "+anchor.getAnalyteMSn().getRt());
    			}
  				count ++;
  			}
  		}
  	}
  	System.out.println(count);
  	System.out.println("data parsed");
  }
  
  @Override
  protected void back() 
  {
  	reset();
  	goBack();
  }
  
  @Override
  protected void next() 
  {
  	boolean noLabelSet = true;
  	boolean deuteratedLabelsPresent = false;
  	for (IsotopeLabelVO labelVO : unambiguousIsotopeLabels_)
  	{
  		if (labelVO.getOmegaPosition() > 0)
  		{
  			noLabelSet = false;
  			if (labelVO.getLabelElements().containsKey(ELEMENT_D))
    		{
    			deuteratedLabelsPresent = true;
    		}
  		}
  	}
  	if (deuteratedLabelsPresent)
  	{
  		try
  		{
  			if (!(getDefaultComponents().getNextPanel() instanceof ExperimentStandardsFileChooserPanel))
  			{
  				getDefaultComponents().addOptionPanelAfterCurrent(new ExperimentStandardsFileChooserPanel(getDefaultComponents(), this));
  			}
  		} catch (Exception ex) {}
  	}
  	else
  	{
  		try
  		{
  			if (getDefaultComponents().getNextPanel() instanceof ExperimentStandardsFileChooserPanel)
  			{
  				getDefaultComponents().removeOptionPanel(getDefaultComponents().getNextPanel());
  			}
  		} catch (Exception ex) {}
  	}
  	if (!noLabelSet)
  	{
  		goNext();
  	}
  	else
  	{
  		new WarningMessage(new JFrame(), "Warning", "You must specify the \u03C9-position for at least one label before continuing!");
  	}
  }
  
  @Override
  protected void updateValue(int row, int value)
	{
  	unambiguousIsotopeLabels_.get(row).setOmegaPosition(value);
	}
  
  public Vector<IsotopeLabelVO> getAssignedIsotopeLabels()
  {
  	Vector<IsotopeLabelVO> assignedIsotopeLabels = new Vector<IsotopeLabelVO>();
  	for (IsotopeLabelVO label : unambiguousIsotopeLabels_)
  	{
  		if (label.getOmegaPosition() > 0)
  		{
  			assignedIsotopeLabels.add(label);
  		}
  	}
  	return assignedIsotopeLabels;
  }
  
  class CalibrationAnchor
  {
  	private String molecularSpecies_;
  	private String filePath_;
  	private LipidomicsMSnSet analyteMSn_;
  	
		public CalibrationAnchor(String molecularSpecies, String filePath,
				LipidomicsMSnSet analyteMSn)
		{
			super();
			this.molecularSpecies_ = molecularSpecies;
			this.filePath_ = filePath;
			this.analyteMSn_ = analyteMSn;
		}
		
		public String getMolecularSpecies()
		{
			return molecularSpecies_;
		}
		
		public String getFilePath()
		{
			return filePath_;
		}
		
		public LipidomicsMSnSet getAnalyteMSn()
		{
			return analyteMSn_;
		}
  }
}
