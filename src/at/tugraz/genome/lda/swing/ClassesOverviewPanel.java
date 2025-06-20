/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2017 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp
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

package at.tugraz.genome.lda.swing;

import java.awt.BorderLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Hashtable;
import java.util.Vector;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JPanel;

import at.tugraz.genome.lda.TooltipTexts;
import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.analysis.ComparativeAnalysis;
import at.tugraz.genome.lda.analysis.SampleLookup;
import at.tugraz.genome.lda.analysis.exception.CalculationNotPossibleException;
import at.tugraz.genome.lda.utils.DoubleCalculator;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.ResultCompGroupVO;
import at.tugraz.genome.lda.vos.ResultCompVO;
import at.tugraz.genome.lda.vos.ResultDisplaySettingsVO;
import at.tugraz.genome.maspectras.utils.Calculator;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class ClassesOverviewPanel extends JPanel implements ActionListener
{

  private static final long serialVersionUID = 7904365893705072363L;
  
  private JCheckBox considerStandard_;
  private JCheckBox considerDilution_;
  private JButton buttonAccept_;
  private JButton buttonSearches_;
  private JPanel centerPanel_;
  
//  private boolean grouped_;
  private Vector<String> expNames_;
  private Hashtable<String,String> fromShortToExpName_;
  private ResultSelectionSettings selectionSettings_;
  private SampleLookup sampleLookup_;
  private Hashtable<String,ResultDisplaySettings> displaySettings_;
  private StandardSelectionOverview standardSelection_;
  private Hashtable<String,Hashtable<String,Integer>> isLookup_;
  private Hashtable<String,Hashtable<String,Integer>> esLookup_;
  private Hashtable<String,Hashtable<String,Hashtable<String,ResultCompVO>>> analysisResults_;
  private Hashtable<String,HeatMapDrawing> selectedIsotopes_;
  private ComparativeAnalysis analysisModule_;
  private ColorChooserDialog colorChooser_;
  
  public ClassesOverviewPanel(Vector<String> expNames,SampleLookup sampleLookup, ComparativeAnalysis analysisModule,
      Hashtable<String,ResultDisplaySettings> displaySettings, Hashtable<String,HeatMapDrawing> selectedIsotopes,
      Hashtable<String,Hashtable<String,Integer>> isLookup, Hashtable<String,Hashtable<String,Integer>> esLookup,
      Hashtable<String,Hashtable<String,Hashtable<String,ResultCompVO>>> analysisResults,ColorChooserDialog colorChooser){
//    grouped_ = grouped;
    expNames_ = expNames;
    sampleLookup_ = sampleLookup;
    displaySettings_ = displaySettings;
    isLookup_ = isLookup;
    esLookup_ = esLookup;
    analysisResults_ = analysisResults;
    selectedIsotopes_ = selectedIsotopes;
    analysisModule_ = analysisModule;
    colorChooser_ = colorChooser;
    init();
  }
  
  private Vector<String> getDisplayNames(){
    Vector<String> names = new Vector<String>();
    fromShortToExpName_ = new Hashtable<String,String>();
    for (String name : this.expNames_){
      String displayName = sampleLookup_.getDisplayName(name);
      names.add(displayName);
      fromShortToExpName_.put(displayName, name);
    }
    return names;
  }
  
  private void init(){
    setLayout(new BorderLayout());
    centerPanel_ = new JPanel();
    add(centerPanel_,BorderLayout.CENTER);
    initBottomPanel();
    selectionSettings_ = new ResultSelectionSettings(null,getDisplayNames(),true);
    selectionSettings_.addActionListener(this);
    standardSelection_ = new StandardSelectionOverview(displaySettings_,this,isLookup_,esLookup_);
  }

  private void initBottomPanel(){
    JPanel bottomPanel = new JPanel();
    bottomPanel.setLayout(new GridBagLayout());
    
    JPanel inputPanel = new JPanel();
    inputPanel.setLayout(new GridBagLayout());
    considerStandard_ = new JCheckBox("consider standard");
    considerStandard_.setSelected(true);
    considerStandard_.setToolTipText(TooltipTexts.OVERVIEW_CONSIDER_STANDARD);
    inputPanel.add(considerStandard_,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    considerDilution_ = new JCheckBox("consider dilution");
    considerDilution_.setSelected(true);
    considerDilution_.setToolTipText(TooltipTexts.OVERVIEW_CONSIDER_DILUTION);
    inputPanel.add(considerDilution_,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0)); 
    bottomPanel.add(inputPanel,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    
    JPanel buttonPanel = new JPanel();
    buttonAccept_ = new JButton("Accept");
    buttonAccept_.addActionListener(this);
    buttonAccept_.setActionCommand("showAcceptOverview");
    buttonAccept_.setToolTipText(TooltipTexts.ACCEPT_GENERAL);
    buttonPanel.add(buttonAccept_,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 2, 0, 2), 0, 0));
    buttonSearches_ = new JButton("Selected");
    buttonSearches_.addActionListener(this);
    buttonSearches_.setActionCommand("showSearches");
    buttonSearches_.setToolTipText(TooltipTexts.OVERVIEW_CHOOSE_SELECTED);
    buttonPanel.add(buttonSearches_,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 2, 0, 2), 0, 0));   
    bottomPanel.add(buttonPanel,new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));    
    add(bottomPanel,BorderLayout.SOUTH);
  }

  public void refreshNames(){
    selectionSettings_.refreshNames(getDisplayNames());
  }
  
  public void actionPerformed(ActionEvent e)
  {
    if (e.getActionCommand().equalsIgnoreCase("showSearches")){
      selectionSettings_.setVisible(true);
    }
    if (e.getActionCommand().equalsIgnoreCase("AcceptDisplaySettings")){
      selectionSettings_.setVisible(false);     
    }
    if (e.getActionCommand().equalsIgnoreCase("showAcceptOverview")){
      if (considerStandard_.isSelected()){
        standardSelection_.refreshInformationSections();
        standardSelection_.setVisible(true);
      } else
        showOverviewBarChart();
//        overviewSelection_.showOverviewPage(grouped_,considerStandard_.isSelected(),considerDilution_.isSelected(),selectionSettings_.getSelected());
    }
    if (e.getActionCommand().equalsIgnoreCase("acceptSelectedStandards")){
      standardSelection_.setVisible(false);
      showOverviewBarChart();
//      overviewSelection_.showOverviewPage(grouped_,considerStandard_.isSelected(),considerDilution_.isSelected(),selectionSettings_.getSelected());
    }
    if (e.getActionCommand().equalsIgnoreCase("cancelSelectedStandards")){
      standardSelection_.setVisible(false);
    }
  }
  
  
  private void showOverviewBarChart(){
    Vector<String> experiments = new Vector<String>();
    for (String displayName : selectionSettings_.getSelected()) experiments.add(fromShortToExpName_.get(displayName));
    Hashtable<String,Hashtable<String,ResultCompVO>> overviewResults = new Hashtable<String,Hashtable<String,ResultCompVO>>();
    Hashtable<String,Hashtable<String,Double>> dilutionFactors = null;
    if (analysisModule_.hasAbsoluteSettings())
      dilutionFactors = analysisModule_.getDilutionFactors();
    Vector<String> classes = new Vector<String>();
    int maxIsotope = Integer.MAX_VALUE;
    boolean isGrouped = false;
    try{
    for (String group : analysisResults_.keySet()){
      for (String molName : analysisResults_.get(group).keySet()){
        for (String expName : analysisResults_.get(group).get(molName).keySet()){
          ResultCompVO compVO = analysisResults_.get(group).get(molName).get(expName);
          if (compVO!=null){
            if (compVO instanceof ResultCompGroupVO){
              isGrouped = true;
              break;
            }else{
              isGrouped = false;
              break;              
            }
          }
        }
      }
    }
    for (String group : analysisResults_.keySet()){
      Hashtable<String,Hashtable<String,ResultCompVO>> resultsOfGroup = analysisResults_.get(group);
      ResultDisplaySettingsVO displaySettings = displaySettings_.get(group).getSettingsVO();
      Hashtable<String,ResultCompVO> overviewOfOneExperiment = new  Hashtable<String,ResultCompVO>();
      if (overviewResults.containsKey(group)) overviewOfOneExperiment = overviewResults.get(group);

      int isotopeGroup = selectedIsotopes_.get(group).getSelectedIsotope()+1;
      if (isotopeGroup>=0 && isotopeGroup<maxIsotope)
        maxIsotope = isotopeGroup;
      boolean isOK = true;
      if (considerStandard_.isSelected()){
//        boolean intern = false;
//        int correctionType = -1;
        if (displaySettings.getESStandMethod()!=ResultCompVO.NO_STANDARD_CORRECTION){
//          correctionType = displaySettings.getESStandMethod();
        }else if (displaySettings.getISStandMethod()!=ResultCompVO.NO_STANDARD_CORRECTION){
//          intern = true;
//          correctionType = displaySettings.getISStandMethod();
        } else isOK = false;
//        if (isOK){
//          double volume = analysisModule_.getStandardAmount(intern, group, correctionType);
//          if (volume>0){
//            refQuantities.put(group, volume);
//          }
//        }
      }

      if (isOK){
        classes.add(group);
        for (String expName : experiments){
        // this experiment cannot be compared if we wish to compare with a standard and there is none selected 
          if (considerStandard_.isSelected()&&displaySettings.getISStandMethod()==ResultCompVO.NO_STANDARD_CORRECTION && displaySettings.getESStandMethod()==ResultCompVO.NO_STANDARD_CORRECTION){
            overviewOfOneExperiment.put(expName, new ResultCompVO(ResultCompVO.CLASS_TYPE,maxIsotope, new Vector<Double>()));
          }else{
            Hashtable<String,Hashtable<Integer,Double>> isoAreas = new Hashtable<String,Hashtable<Integer,Double>>();
            for (int isoNr=0; isoNr!=(isotopeGroup+1); isoNr++){
              for (String molName : resultsOfGroup.keySet()){
                Hashtable<String,ResultCompVO> resultsOfMol = resultsOfGroup.get(molName);
                if (resultsOfMol.containsKey(expName)){
                  ResultCompVO compVO = resultsOfMol.get(expName);
                  if (isGrouped){
                    Hashtable<String,ResultCompVO> groupVOs = ((ResultCompGroupVO)compVO).getGroupingPartners();
                    for (String origExp : groupVOs.keySet()){
                      Hashtable<Integer,Double> nrHash = new Hashtable<Integer,Double>();
                      if (isoAreas.containsKey(origExp)) nrHash = isoAreas.get(origExp);
                      double groupArea = 0;
                      if (nrHash.containsKey(isoNr)) groupArea = nrHash.get(isoNr);
                      groupArea += getCorrespondingArea(groupVOs.get(origExp), displaySettings,isoNr,dilutionFactors,group,origExp);
//                      System.out.println(group+"/"+expName+"/"+"/"+molName+"/"+isoNr+"/"+groupArea);
                      nrHash.put(isoNr, groupArea);
                      isoAreas.put(origExp, nrHash);
                    }
                  }else{
                    Hashtable<Integer,Double> nrHash = new Hashtable<Integer,Double>();
                    if (isoAreas.containsKey(expName)) nrHash = isoAreas.get(expName);
                    double groupArea = 0;
                    if (nrHash.containsKey(isoNr)) groupArea = nrHash.get(isoNr);
                    groupArea += getCorrespondingArea(compVO, displaySettings,isoNr,dilutionFactors,group,expName);
                    nrHash.put(isoNr, groupArea);
                    isoAreas.put(expName, nrHash);
                  }
                }
              }
            }
            ResultCompVO compVO;
            if (isGrouped){
              Hashtable<String,ResultCompVO> groupMembers = new Hashtable<String,ResultCompVO>();
              for (String origExpName : isoAreas.keySet()){
                Hashtable<Integer,Double> values = isoAreas.get(origExpName);
                Vector<Double> areas = new Vector<Double>();
                for (int isoNr=0; isoNr!=(isotopeGroup+1); isoNr++) areas.add(values.get(isoNr));
                groupMembers.put(origExpName, new ResultCompVO(ResultCompVO.CLASS_TYPE,maxIsotope,areas));
              }
              compVO = new ResultCompGroupVO(groupMembers);
            }else{
              Vector<Double> areas = new Vector<Double>();
              Hashtable<Integer,Double> values = isoAreas.values().iterator().next();
              for (int isoNr=0; isoNr!=(isotopeGroup+1); isoNr++) areas.add(values.get(isoNr));
              compVO = new ResultCompVO(ResultCompVO.CLASS_TYPE,maxIsotope,areas);              
            }
            //overviewOfOneExperiment.put(sampleLookup_.getDisplayName(expName), compVO);
            overviewOfOneExperiment.put(expName, compVO);
          }
//          overviewResults.put(expName, overviewOfOneExperiment);
        }
        overviewResults.put(group, overviewOfOneExperiment);
      }    
    }
    if (maxIsotope==Integer.MAX_VALUE)
      maxIsotope = 1;
    Hashtable<String,Vector<Double>> totalAmounts = new Hashtable<String,Vector<Double>>();
    Hashtable<String,Vector<Double>> totalAmountsOfExp = new Hashtable<String,Vector<Double>>();
    Hashtable<String,Vector<Double>> totalSds = new Hashtable<String,Vector<Double>>();
    Vector<Double> allValues = new Vector<Double>();
    //for (String expName : selectionSettings_.getSelected()){
    for (String expName : experiments){  
      Vector<Double> totalAreas = new Vector<Double>();
      Vector<Double> sds = new Vector<Double>();
      Hashtable<String,Vector<Double>> singleAreasIso = new Hashtable<String,Vector<Double>>();
      for (int isoNr=0; isoNr!=(maxIsotope+1); isoNr++){
        double area = 0;
        Hashtable<String,Double> singleExpAreas = new Hashtable<String,Double>();
        Vector<Double> stdvs = new Vector<Double>();
        for (String group : analysisResults_.keySet()){
          if (overviewResults.containsKey(group) && overviewResults.get(group).containsKey(expName)){
            ResultCompVO compVO = overviewResults.get(group).get(expName);
            double value = compVO.getOriginalArea(compVO.getAvailableIsotopeNr(isoNr));
            area += value;
            if (compVO instanceof ResultCompGroupVO){
              ResultCompGroupVO groupVO = (ResultCompGroupVO)compVO;
              stdvs.add(groupVO.getOriginalAreaSD(compVO.getAvailableIsotopeNr(isoNr)));
              for (String origExpName : groupVO.getGroupingPartners().keySet()){
                double singleExpArea = 0;
                if (singleExpAreas.containsKey(origExpName)) singleExpArea = singleExpAreas.get(origExpName);
                singleExpArea += groupVO.getGroupingPartners().get(origExpName).getOriginalArea(compVO.getAvailableIsotopeNr(isoNr));
                singleExpAreas.put(origExpName, singleExpArea);
              }
            }
            if (value>0 && isoNr==maxIsotope)
              allValues.add(value);
          }
        }
        totalAreas.add(area);
        if (isGrouped){
          sds.add(Calculator.calculateSumStdevErrorPropagated(stdvs));
          for (String origExp : singleExpAreas.keySet()){
            Vector<Double> areas = new Vector<Double>();
            if (singleAreasIso.containsKey(origExp)) areas = singleAreasIso.get(origExp);
            areas.add(singleExpAreas.get(origExp));
            singleAreasIso.put(origExp, areas);
          }
        }
      }
      totalAmounts.put(expName, totalAreas);
      if (isGrouped){
        totalSds.put(expName, sds);
        for (String origExp : singleAreasIso.keySet()){
          totalAmountsOfExp.put(origExp,singleAreasIso.get(origExp));
        }
      }
    }
    String preferredUnit = "";
    String unit = "area [AU]";
    if (considerStandard_.isSelected() && analysisModule_.hasAbsoluteSettings()){
      double median = DoubleCalculator.median(allValues);
      preferredUnit = StaticUtils.extractPreferredUnit(median);
      unit = StaticUtils.getCorrespondingUnit(new ResultDisplaySettingsVO("amount sample-volume",ResultCompVO.NO_STANDARD_CORRECTION,ResultCompVO.NO_STANDARD_CORRECTION,false,false), preferredUnit,true);
    }
    for (String group : overviewResults.keySet()){
      for (String exp : overviewResults.get(group).keySet()){        
        overviewResults.get(group).get(exp).setSumValueForPercentage(totalAmounts.get(exp));
        if (isGrouped){
          ResultCompGroupVO compVO = (ResultCompGroupVO)overviewResults.get(group).get(exp);
          compVO.setSumPercentualMeans(totalAmounts.get(exp));
          compVO.setSumPercentualSds(totalSds.get(exp));
          Hashtable<String,ResultCompVO> hashes = ((ResultCompGroupVO)overviewResults.get(group).get(exp)).getGroupingPartners();
          for (String origExp : hashes.keySet()){
            ResultCompVO groupVO = hashes.get(origExp);
            groupVO.setSumValueForPercentage(totalAmountsOfExp.get(origExp));
          }
        }
      }
    }
    
    ResultDisplaySettingsVO settingVO = new ResultDisplaySettingsVO(ResultDisplaySettingsVO.REL_VALUE, ResultCompVO.NO_STANDARD_CORRECTION, ResultCompVO.NO_STANDARD_CORRECTION,
        false, false);
    centerPanel_.removeAll();
    BarChartPainter painter = new BarChartPainter(BarChartPainter.TYPE_CLASSES,"Overview",classes,overviewResults,experiments,sampleLookup_,true,isGrouped,maxIsotope,false,isGrouped,settingVO, preferredUnit,unit,
        null,null,null,colorChooser_);
    this.remove(centerPanel_);
    this.centerPanel_ = painter;
    this.add(centerPanel_);
    this.invalidate();
    this.updateUI();
    }catch(CalculationNotPossibleException cnp){
      new WarningMessage(new JFrame(), "Error", cnp.getMessage());
    }
//    Hashtable<String,String> fromShortNameToExperiment = new Hashtable<String,String>();
//    for (String selectionSettings_.getSelected()
  }
  
  private double getCorrespondingArea(ResultCompVO compVO, ResultDisplaySettingsVO displaySettings, int isoNr, Hashtable<String,Hashtable<String,Double>> dilutionFactors, String group, String expName) throws CalculationNotPossibleException{
    double area = 0;
    if (compVO.getType()==ResultCompVO.ANALYTE_TYPE){
      if (considerStandard_.isSelected()){
//        ResultDisplaySettingsVO settingVO = new ResultDisplaySettingsVO(null,displaySettings.getISStandMethod(), displaySettings.getESStandMethod(), considerDilution_.isSelected(),false);
//        area += compVO.getRelativeToMedian(compVO.getAvailableIsotopeNr(isoNr), settingVO)*refQuantities.get(group);
        double amount = compVO.getAmountInProbeVolume(compVO.getAvailableIsotopeNr(isoNr), displaySettings.getISStandMethod(), displaySettings.getESStandMethod());
        if (!considerDilution_.isSelected() && dilutionFactors!=null && dilutionFactors.containsKey(group) && dilutionFactors.get(group).containsKey(expName)){
          amount =  amount/dilutionFactors.get(group).get(expName);
        }  
        area += amount;
      }else{
        area += compVO.getStandardizedArea(compVO.getAvailableIsotopeNr(isoNr), ResultCompVO.NO_STANDARD_CORRECTION, ResultCompVO.NO_STANDARD_CORRECTION, considerDilution_.isSelected());
      }
    }  
    return area;
  }
  
  public void cleanup()
  {
  	if (selectionSettings_!=null) 
  	{
  		selectionSettings_.removeActionListener(this);
  		selectionSettings_.dispose();
    	selectionSettings_ = null;
  	}
  	if (standardSelection_!=null) 
  	{
  		standardSelection_.cleanup();
    	standardSelection_ = null;
  	}	
  }
}
