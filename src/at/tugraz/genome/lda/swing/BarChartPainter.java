/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2017 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger 
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
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusEvent;
import java.awt.event.ItemEvent;
import java.awt.event.KeyEvent;
import java.awt.event.FocusListener;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.lang.reflect.Method;
import java.util.Hashtable;
import java.util.Vector;

import javax.imageio.ImageIO;
import javax.swing.ButtonGroup;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JTextField;
import javax.swing.InputVerifier;
import javax.swing.filechooser.FileNameExtensionFilter;

import org.apache.batik.dom.GenericDOMImplementation;
import org.apache.batik.svggen.SVGGraphics2D;
import org.w3c.dom.DOMImplementation;
import org.w3c.dom.Document;

import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.TooltipTexts;
import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.analysis.SampleLookup;
import at.tugraz.genome.lda.analysis.exception.CalculationNotPossibleException;
import at.tugraz.genome.lda.exception.ExcelInputFileException;
import at.tugraz.genome.lda.exception.ExportException;
import at.tugraz.genome.lda.export.ExcelAndTextExporter;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.ExportOptionsVO;
import at.tugraz.genome.lda.vos.ResultCompGroupVO;
import at.tugraz.genome.lda.vos.ResultCompVO;
import at.tugraz.genome.lda.vos.ResultDisplaySettingsVO;
import at.tugraz.genome.maspectras.parser.exceptions.SpectrummillParserException;
import at.tugraz.genome.maspectras.utils.Calculator;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class BarChartPainter extends JPanel implements ActionListener
{
  
  private static final long serialVersionUID = -8728298412578227539L;
  
  public final static int TYPE_MOLECULE = 0;
  public final static int TYPE_EXPERIMENT = 1;
  public final static int TYPE_CLASSES = 2;

  private Vector<String> moleculeNames_;
  private JComboBox<String> quantType_;
  private JPanel barChartPanel_;
  
  private ButtonGroup barDisplayGroup_;
  private JRadioButton singleSided_;
  private JRadioButton doubleSided_;
  
  private ButtonGroup barScaleGroup_;
  private JRadioButton logarithmicScale_;
  private JRadioButton normalScale_;
  
  private JComboBox<String> magnitudeChooser_;
  
  
  private BarChart barChart_;
  
  private JRadioButton standardDeviation_;
  private JTextField multStandardDeviation_;
  private JRadioButton standardErrorMean_;
  protected JCheckBox exportDeviation_;
  private JRadioButton columnAnalyte_;
  private JRadioButton columnExperiment_;
  protected JCheckBox exportRT_;
  protected JCheckBox exportRTDev_;

  protected JComboBox<String> maxIsotopes_;
  
  private Hashtable<String,Hashtable<String,ResultCompVO>> analysisResultsHash_;
  
  private boolean repaintingBarChart_;
  private Vector<String> originalValueNames_;
  private boolean showSDFields_;
  private int maximumIsotopes_;
  private ResultDisplaySettingsVO settingVO_;
  private String preferredUnit_;
  private String unit_;
  private JFileChooser exportFileChooser_;
  
  private int type_;
  private String groupName_;
  
  private Hashtable<String,Hashtable<String,Double>> valuesToPaint_ = new Hashtable<String,Hashtable<String,Double>>();
  private Hashtable<String,Hashtable<String,Double>> sdsToPaint_ = new Hashtable<String,Hashtable<String,Double>>();  
  /** the results from the original experiments - the mean and SD values are calculated in the new exporter; first key: analyte name; second key experiment name*/
  private Hashtable<String,Hashtable<String,Double>> singleResultValues_ = new Hashtable<String,Hashtable<String,Double>>();
  private Hashtable<String,Integer> corrTypeISLookup_;
  private Hashtable<String,Integer> corrTypeESLookup_;
  Vector<String> modifications_;
  private String yAxisText_;
  ColorChooserDialog colorChooser_;
  private boolean bigConstructor_ = false;
  private int colorType_;
  private SampleLookup sampleLookup_;
  /** should for this instance the sample lookupup be used*/
  private boolean useSampleLookup_;
  /** true when the peaks were grouped by a certain retention time*/
  private boolean rtGrouped_;
  /** true when sample groups are exported*/
  private boolean isGroupedView_;
  
  /**
   * constructor for a bar chart painter when a single analyte/experiment/sample group is clicked
   * @param type the export type (see TYPE_ ... definitions in the variables)
   * @param groupName the analyte class
   * @param moleculeName the name of the analyte/experiment/sample 
   * @param analysisResults the corresponding results
   * @param sortedValueNames the names of x-axis values
   * @param lookup a lookup to the original sample names
   * @param useSampleLookup true when the lookup for the original sample names shall be used
   * @param showSDFields should error bars be shown
   * @param maxIsotopes the maximum isotope that will be used for the export
   * @param rtGrouped true when the peaks were grouped by a certain retention time
   * @param isGroupedView true when sample groups are exported
   * @param settingVO value object specifying the type of displayed values
   * @param preferredUnit the preferred unit magnifier, e.g p for pico
   * @param unit a description of the physical unit
   * @param corrTypeISLookup the internal standard correction type
   * @param corrTypeESLookup the external standard correction type
   * @param modifications the modifications for this analyte class
   * @param colorChooser the color chooser dialog
   */
  public BarChartPainter(int type, String groupName, String moleculeName,Hashtable<String,ResultCompVO> analysisResults,Vector<String> sortedValueNames, SampleLookup lookup,
      boolean useSampleLookup, boolean showSDFields, int maxIsotopes, boolean rtGrouped, boolean isGroupedView, ResultDisplaySettingsVO settingVO, String preferredUnit, String unit,
      Hashtable<String,Integer> corrTypeISLookup, Hashtable<String,Integer> corrTypeESLookup, Vector<String> modifications, ColorChooserDialog colorChooser){
    this(type, groupName, lookup, useSampleLookup, sortedValueNames, showSDFields, maxIsotopes, rtGrouped, isGroupedView, settingVO, preferredUnit,
        unit, corrTypeISLookup, corrTypeESLookup, modifications, colorChooser);
    moleculeNames_.add(moleculeName);
    analysisResultsHash_.put(moleculeName,analysisResults);
    bigConstructor_ = false;
    init();
  }
  
  /**
   * constructor for a bar chart painter when there are several molecules selected
   * @param type the export type (see TYPE_ ... definitions in the variables)
   * @param groupName the analyte class
   * @param moleculeNames the names of the analytes
   * @param analysisResults the corresponding results
   * @param sortedValueNames the names of x-axis values
   * @param lookup a lookup to the original sample names
   * @param useSampleLookup true when the lookup for the original sample names shall be used
   * @param showSDFields should error bars be shown
   * @param maxIsotopes the maximum isotope that will be used for the export
   * @param rtGrouped true when the peaks were grouped by a certain retention time
   * @param isGroupedView true when sample groups are exported
   * @param settingVO value object specifying the type of displayed values
   * @param preferredUnit the preferred unit magnifier, e.g p for pico
   * @param unit a description of the physical unit
   * @param corrTypeISLookup the internal standard correction type
   * @param corrTypeESLookup the external standard correction type
   * @param modifications the modifications for this analyte class
   * @param colorChooser the color chooser dialog
   */
  public BarChartPainter(int type, String groupName, Vector<String> moleculeNames,Hashtable<String,Hashtable<String,ResultCompVO>> analysisResults,Vector<String> sortedValueNames,
      SampleLookup lookup, boolean useSampleLookup, boolean showSDFields, int maxIsotopes, boolean rtGrouped, boolean isGroupedView, ResultDisplaySettingsVO settingVO,
      String preferredUnit, String unit,Hashtable<String,Integer> corrTypeISLookup, Hashtable<String,Integer> corrTypeESLookup, Vector<String> modifications,
      ColorChooserDialog colorChooser){
    this(type, groupName, lookup, useSampleLookup, sortedValueNames, showSDFields, maxIsotopes, rtGrouped, isGroupedView, settingVO, preferredUnit,
        unit,corrTypeISLookup, corrTypeESLookup, modifications, colorChooser);
    analysisResultsHash_ = analysisResults;
    moleculeNames_ = moleculeNames;
    bigConstructor_ = true;
    init();
  }
  
  /**
   * constructor that have the other constructors in common
   * @param type the export type (see TYPE_ ... definitions in the variables)
   * @param groupName the analyte class
   * @param lookup a lookup to the original sample names
   * @param useSampleLookup true when the lookup for the original sample names shall be used
   * @param sortedValueNames the names of x-axis values
   * @param showSDFields should error bars be shown
   * @param maxIsotopes the maximum isotope that will be used for the export
   * @param rtGrouped true when the peaks were grouped by a certain retention time
   * @param isGroupedView true when sample groups are exported 
   * @param settingVO value object specifying the type of displayed values
   * @param preferredUnit the preferred unit magnifier, e.g p for pico
   * @param unit a description of the physical unit
   * @param corrTypeISLookup the internal standard correction type
   * @param corrTypeESLookup the external standard correction type
   * @param modifications the modifications for this analyte class
   * @param colorChooser the color chooser dialog
   */
  private BarChartPainter(int type, String groupName, SampleLookup lookup, boolean useSampleLookup, Vector<String> sortedValueNames, boolean showSDFields,
      int maxIsotopes, boolean rtGrouped, boolean isGroupedView, ResultDisplaySettingsVO settingVO, String preferredUnit, String unit,
      Hashtable<String,Integer> corrTypeISLookup,  Hashtable<String,Integer> corrTypeESLookup, Vector<String> modifications,
      ColorChooserDialog colorChooser){
    analysisResultsHash_ = new Hashtable<String,Hashtable<String,ResultCompVO>>();
    sampleLookup_ = lookup;
    useSampleLookup_ = useSampleLookup;
    moleculeNames_ = new Vector<String>();
    this.type_ = type;
    this.groupName_ = groupName;
    originalValueNames_ = sortedValueNames;
    maximumIsotopes_ = maxIsotopes;
    rtGrouped_ = rtGrouped;
    isGroupedView_ = isGroupedView;
    showSDFields_ = showSDFields;
    settingVO_ = settingVO;
    preferredUnit_ = preferredUnit;
    unit_ = unit;
    modifications_ = modifications;
    if (modifications_==null) modifications_=new Vector<String>();
    this.setLayout(new BorderLayout());
    barChartPanel_ = new JPanel();
    barChartPanel_.setLayout(new BorderLayout());
    this.add(barChartPanel_,BorderLayout.CENTER);
    corrTypeISLookup_ = corrTypeISLookup;
    corrTypeESLookup_ = corrTypeESLookup;
    colorChooser_ = colorChooser;
  }
  
  private void init(){
    colorType_ = ColorChooserDialog.DEFAULT_TYPE;
    if (bigConstructor_){
      colorType_ = ColorChooserDialog.EXPERIMENT_TYPE;
      if (showSDFields_)
        colorType_ = ColorChooserDialog.GROUP_TYPE;
    }
    this.createTopMenu();
    ExportPanel expPanel = new ExportPanel(Color.WHITE,Color.BLACK,this);
    this.add(expPanel,BorderLayout.SOUTH);
    repaintingBarChart_ = false;
    exportFileChooser_ = new JFileChooser();
    exportFileChooser_.setPreferredSize(new Dimension(600,500));
    
  }
  
  private void createTopMenu(){
    JPanel menuPanel = new JPanel();
    menuPanel.setLayout(new BorderLayout());
    this.add(menuPanel,BorderLayout.NORTH);
    Font titleFont=new Font("Helvetica",Font.PLAIN,24);
    String titleName = "";
    if (this.moleculeNames_.size()<2)
      titleName = moleculeNames_.get(0);
    JLabel title = new JLabel(titleName);
    title.setFont(titleFont);
    title.setHorizontalAlignment(JLabel.CENTER);
    menuPanel.add(title,BorderLayout.NORTH);

    maxIsotopes_ = new JComboBox<String>();
    for (int i=0;i!=maximumIsotopes_;i++){
      maxIsotopes_.addItem(String.valueOf(i));
    }
    maxIsotopes_.setSelectedItem(String.valueOf(maximumIsotopes_-1));

    
    quantType_ = new JComboBox<String>();
    quantType_.addItem("area absolute");
    quantType_.addItem("percentual values");
    if (type_ != TYPE_CLASSES)
      quantType_.addItem("area relative to molecule");
    if (settingVO_.getISStandMethod() != ResultCompVO.NO_STANDARD_CORRECTION || settingVO_.getESStandMethod() != ResultCompVO.NO_STANDARD_CORRECTION)
      quantType_.addItem("area relative to standard");
    if (!settingVO_.getType().equalsIgnoreCase("relative value"))
      quantType_.addItem("absolute quantity");
    repaintingBarChart_ = true;
    
    if (type_ != TYPE_CLASSES)
      quantType_.setSelectedItem("area relative to molecule");
    else
      quantType_.setSelectedItem("percentual values");
    repaintingBarChart_ = false;
    SelectionItemListener quantTypeListener = new SelectionItemListener("ChangeChart");
    quantType_.addItemListener(quantTypeListener);
    quantType_.setToolTipText(TooltipTexts.BARCHART_VALUE_TYPE);
    
    magnitudeChooser_ = new JComboBox<String>();
    addEmptyItemToChooser();
    magnitudeChooser_.setEnabled(false);
    SelectionItemListener magnListener = new SelectionItemListener("ChangeMagnitude");
    magnitudeChooser_.addItemListener(magnListener);
    magnitudeChooser_.setToolTipText(TooltipTexts.BARCHART_MAGNITUDE_CHOOSER);
    
    barDisplayGroup_ = new ButtonGroup();
    singleSided_ = new JRadioButton("single-sided");
    singleSided_.addItemListener(new SelectionItemListener("ChangeSided"));
    singleSided_.setToolTipText(TooltipTexts.BARCHART_SINGLE_SIDED);
    barDisplayGroup_.add(singleSided_); 
    doubleSided_ = new JRadioButton("double-sided");
    doubleSided_.addItemListener(new SelectionItemListener("ChangeSided"));
    doubleSided_.setToolTipText(TooltipTexts.BARCHART_DOUBLE_SIDED);
    barDisplayGroup_.add(doubleSided_); 
    
    barScaleGroup_ = new ButtonGroup();
    logarithmicScale_ = new JRadioButton("logarithmic");
    logarithmicScale_.addItemListener(new SelectionItemListener("ChangeLog"));
    logarithmicScale_.setToolTipText(TooltipTexts.BARCHART_LOGARITHMIC);
    barScaleGroup_.add(logarithmicScale_);
    normalScale_ = new JRadioButton("linear");
    normalScale_.setSelected(true);
    normalScale_.addItemListener(new SelectionItemListener("ChangeLog"));
    normalScale_.setToolTipText(TooltipTexts.BARCHART_LINEAR);
    barScaleGroup_.add(normalScale_); 

    SelectionItemListener isotopeListener = new SelectionItemListener("ChangeIsotope");
    maxIsotopes_.addItemListener(isotopeListener);
    maxIsotopes_.setToolTipText(TooltipTexts.HEATMAP_ISOTOPES);
    
    JPanel menuItemsBar = new JPanel();
    JLabel typeLabel = new JLabel("quant-type: ");
    typeLabel.setToolTipText(TooltipTexts.BARCHART_VALUE_TYPE);
    menuItemsBar.add(typeLabel);
    
    menuItemsBar.add(quantType_);
    menuItemsBar.add(magnitudeChooser_);
    menuItemsBar.add(singleSided_);    
    menuItemsBar.add(doubleSided_);
    menuItemsBar.add(logarithmicScale_);    
    menuItemsBar.add(normalScale_);
    menuItemsBar.add(maxIsotopes_);
    JLabel isotpesLabel = new JLabel("isotopes");
    isotpesLabel.setToolTipText(TooltipTexts.HEATMAP_ISOTOPES);
    menuItemsBar.add(isotpesLabel);

    if (settingVO_.getType().equalsIgnoreCase("relative value"))
      quantType_.setSelectedItem("area absolute");
    else
      quantType_.setSelectedItem("absolute quantity");
    menuPanel.add(menuItemsBar,BorderLayout.CENTER);
    JPanel expAndSD = new JPanel();
    expAndSD.setLayout(new BorderLayout());
    menuPanel.add(expAndSD,BorderLayout.SOUTH);
    if (showSDFields_){  
      JPanel sdItemsBar = new JPanel();
      expAndSD.add(sdItemsBar,BorderLayout.NORTH);
      ButtonGroup deviationGroup = new ButtonGroup();
      standardDeviation_ = new JRadioButton("standard deviation");
      standardDeviation_.setSelected(true);
      standardDeviation_.addItemListener(new SelectionItemListener("ChangeSD"));
      standardDeviation_.setToolTipText(TooltipTexts.BARCHART_STANDARD_DEVIATION);
      deviationGroup.add(standardDeviation_);
      sdItemsBar.add(standardDeviation_);
      multStandardDeviation_ = new JTextField(3);
      multStandardDeviation_.setText("1.0");
      multStandardDeviation_.addKeyListener(new ReturnKeyListener());
      multStandardDeviation_.setInputVerifier(new DoubleVerifier());
      multStandardDeviation_.addFocusListener(new SDFocusListener());
      multStandardDeviation_.setToolTipText(TooltipTexts.BARCHART_STANDARD_DEVIATION_WHICH);
      sdItemsBar.add(multStandardDeviation_);
      standardErrorMean_ = new JRadioButton("standard error mean");
      standardErrorMean_.addItemListener(new SelectionItemListener("ChangeSD"));
      standardErrorMean_.setToolTipText(TooltipTexts.BARCHART_STANDARD_ERROR);
      deviationGroup.add(standardErrorMean_);
      sdItemsBar.add(standardErrorMean_);
    }
    JPanel expItemsBar = new JPanel();
    expAndSD.add(expItemsBar,BorderLayout.SOUTH);
    if (showSDFields_){  
      exportDeviation_ = new JCheckBox("export deviation value"); 
      exportDeviation_.setSelected(false);
      exportDeviation_.setToolTipText(TooltipTexts.HEATMAP_EXPORT_DEVIATION);
      expItemsBar.add(exportDeviation_);
    }
    ButtonGroup columnGroup = new ButtonGroup();
    columnAnalyte_ = new JRadioButton("analyte in column");
    columnAnalyte_.setSelected(true);
    columnAnalyte_.setToolTipText(TooltipTexts.HEATMAP_EXPORT_COLUMN_ANALYTE);
    columnGroup.add(columnAnalyte_);
    expItemsBar.add(columnAnalyte_);    
    columnExperiment_ = new JRadioButton("exp in column");
    columnExperiment_.setSelected(false);
    columnExperiment_.setToolTipText(TooltipTexts.HEATMAP_EXPORT_COLUMN_EXPERIMENT);
    columnGroup.add(columnExperiment_);
    expItemsBar.add(columnExperiment_);
    if (type_!=TYPE_CLASSES){
      exportRT_ = new JCheckBox("export retention-time"); 
      exportRT_.setActionCommand(ExportSettingsPanel.CHANGE_RT_SELECTION_STATUS);
      exportRT_.addActionListener(this);
      exportRT_.setSelected(false);
      exportRT_.setToolTipText(TooltipTexts.HEATMAP_EXPORT_RT);
      expItemsBar.add(exportRT_);
      if (showSDFields_){  
        exportRTDev_ = new JCheckBox("export RT-stdev"); 
        exportRTDev_.addActionListener(this);
        exportRTDev_.setSelected(false);
        exportRTDev_.setEnabled(false);
        exportRTDev_.setToolTipText(TooltipTexts.HEATMAP_EXPORT_RT_SD);
        expItemsBar.add(exportRTDev_); 
      }  
    }
  }
  
  
  @SuppressWarnings("unchecked")
  public void createPaintingArea(){
    barChartPanel_.removeAll();
    valuesToPaint_ = new Hashtable<String,Hashtable<String,Double>>();
    sdsToPaint_ = new Hashtable<String,Hashtable<String,Double>>();
    double lowestValue = 0;
    double highestValue = 0;
    int maxIsotope = Integer.parseInt((String)maxIsotopes_.getSelectedItem());
    yAxisText_ = null;
    String subTitle = "";
    String standardText = "";
    if (settingVO_.getISStandMethod()!=ResultCompVO.NO_STANDARD_CORRECTION || settingVO_.getESStandMethod()!=ResultCompVO.NO_STANDARD_CORRECTION){
      standardText = "standardized:   ";
      if (settingVO_.getISStandMethod()!=ResultCompVO.NO_STANDARD_CORRECTION){
        standardText +="IS=";
        if (settingVO_.getISStandMethod()==ResultCompVO.INTERNAL_STANDARD_TYPE)
          standardText +="most reliable standard";
        else if (settingVO_.getISStandMethod()==ResultCompVO.EXTERNAL_STANDARD_TYPE)
          standardText += "median method";
        else
          standardText += getStandardName(this.corrTypeISLookup_,settingVO_.getISStandMethod());
      }
      if (settingVO_.getISStandMethod()!=ResultCompVO.NO_STANDARD_CORRECTION && settingVO_.getESStandMethod()!=ResultCompVO.NO_STANDARD_CORRECTION)
        standardText += "   ";
      if (settingVO_.getESStandMethod()!=ResultCompVO.NO_STANDARD_CORRECTION){
        standardText +="ES=";
        if (settingVO_.getESStandMethod()==ResultCompVO.INTERNAL_STANDARD_TYPE)
          standardText +="most reliable standard";
        else if (settingVO_.getESStandMethod()==ResultCompVO.EXTERNAL_STANDARD_TYPE)
          standardText += "median method";
        else
          standardText += getStandardName(this.corrTypeESLookup_,settingVO_.getESStandMethod());
      }
      
    }
    if (((String)quantType_.getSelectedItem()).equalsIgnoreCase("area absolute")){
      ResultDisplaySettingsVO standardizedSet = new ResultDisplaySettingsVO(settingVO_);
      standardizedSet.setType("relative value");
      boolean useAbsolute = false;
      if (type_ == TYPE_CLASSES){
        useAbsolute = true;
      }  
      @SuppressWarnings("rawtypes")
      Vector result = this.extractLowestHighestValue("getArea","getAreaSD","getAreaSE",maxIsotope,standardizedSet,useAbsolute,preferredUnit_);
      lowestValue = (Double)result.get(0);
      highestValue = (Double)result.get(1);
      valuesToPaint_ = (Hashtable<String,Hashtable<String,Double>>)result.get(2);
      sdsToPaint_ = (Hashtable<String,Hashtable<String,Double>>)result.get(3);
      singleResultValues_ = (Hashtable<String,Hashtable<String,Double>>)result.get(4);
      yAxisText_ = "area [AU]";
      if (type_ == TYPE_CLASSES){
        yAxisText_ = unit_;
      }else{
        subTitle = "arbitrary area value";
      }
    }
    if (((String)quantType_.getSelectedItem()).equalsIgnoreCase("percentual values")){
      try{
        ResultDisplaySettingsVO standardizedSet = new ResultDisplaySettingsVO(settingVO_);
        subTitle = "relative to total amount of ";
        if ((type_ == TYPE_MOLECULE && bigConstructor_) || type_ == TYPE_EXPERIMENT)
          subTitle += "selected molecules";
        else if (type_ == TYPE_MOLECULE )
          subTitle += "selected samples";
        else if (type_ == TYPE_CLASSES)
          subTitle += "classes";
        if (type_ != TYPE_CLASSES){
          calculatePercentualValues(maxIsotope);
          standardizedSet.setPercent(true);
        }else
          standardizedSet.setType("percentual value");
        standardizedSet.setDivisorMagnitude("");
        @SuppressWarnings("rawtypes")
		Vector result = this.extractLowestHighestValue("getArea","getAreaSD","getAreaSE",maxIsotope,standardizedSet,true,"%");
        lowestValue = (Double)result.get(0);
        highestValue = (Double)result.get(1);
        valuesToPaint_ = (Hashtable<String,Hashtable<String,Double>>)result.get(2);
        sdsToPaint_ = (Hashtable<String,Hashtable<String,Double>>)result.get(3);
        singleResultValues_ = (Hashtable<String,Hashtable<String,Double>>)result.get(4);
        yAxisText_ = "percent [%]";
        
        
      }catch (CalculationNotPossibleException cnp){
        new WarningMessage(new JFrame(),"Error",cnp.getMessage());
        return;
      }
      
    }
    if (((String)quantType_.getSelectedItem()).equalsIgnoreCase("area relative to molecule")){
      @SuppressWarnings("rawtypes")
      Vector result = this.extractLowestHighestValue("getRelativeValue","getRelativeValueSD","getRelativeValueSE",maxIsotope,settingVO_,false,preferredUnit_);
      lowestValue = (Double)result.get(0);
      highestValue = (Double)result.get(1);
      valuesToPaint_ = (Hashtable<String,Hashtable<String,Double>>)result.get(2);
      sdsToPaint_ = (Hashtable<String,Hashtable<String,Double>>)result.get(3);
      singleResultValues_ = (Hashtable<String,Hashtable<String,Double>>)result.get(4);
      yAxisText_ = "rel. median";
      subTitle = "relative to median";
    }
    if (((String)quantType_.getSelectedItem()).equalsIgnoreCase("area relative to standard")){
      @SuppressWarnings("rawtypes")
      Vector result = this.extractLowestHighestValue("getRelativeToMedian","getRelativeToMedianSD","getRelativeToMedianSE",maxIsotope,settingVO_,false,preferredUnit_);
      lowestValue = (Double)result.get(0);
      highestValue = (Double)result.get(1);
      valuesToPaint_ = (Hashtable<String,Hashtable<String,Double>>)result.get(2);
      sdsToPaint_ = (Hashtable<String,Hashtable<String,Double>>)result.get(3);
      singleResultValues_ = (Hashtable<String,Hashtable<String,Double>>)result.get(4);
      yAxisText_ = "rel. standard";
      subTitle = "relative to standard";
    }
    if (lowestValue==1&&highestValue==1){
      highestValue = 10d;
      lowestValue = 0.1d;
    }  
    String sdText = "";
    if (showSDFields_){
      if (this.standardDeviation_!=null){
        if (standardDeviation_.isSelected()){
          sdText = "+/- "+multStandardDeviation_.getText()+"*SD";
        }else{
          sdText = "+/- SE";
        }
      }else{
        sdText = "+/- 1.0*SD";
      }
    }
    if (((String)quantType_.getSelectedItem()).equalsIgnoreCase("absolute quantity")){
      @SuppressWarnings("rawtypes")
      Vector result = this.extractLowestHighestValue("getArea","getAreaSD","getAreaSE",maxIsotope,settingVO_,true,(String)magnitudeChooser_.getSelectedItem());
      lowestValue = (Double)result.get(0);
      highestValue = (Double)result.get(1);
      valuesToPaint_ = (Hashtable<String,Hashtable<String,Double>>)result.get(2);
      sdsToPaint_ = (Hashtable<String,Hashtable<String,Double>>)result.get(3);
      singleResultValues_ = (Hashtable<String,Hashtable<String,Double>>)result.get(4);
      yAxisText_ = unit_.substring(0,unit_.indexOf("[")+1)+(String)magnitudeChooser_.getSelectedItem()+unit_.substring(unit_.indexOf("[")+1+preferredUnit_.length());
      
      if (settingVO_.getType().equalsIgnoreCase("relative to base peak")){
        subTitle = "relative to highest peak of the "+groupName_+" class";
      }else if (settingVO_.getType().equalsIgnoreCase("relative to measured class amount")){
        subTitle = "relative to the total amount of "+groupName_+" class";
      }else if (settingVO_.getType().equalsIgnoreCase("relative to highest total peak")){
        subTitle = "relative to the highest found peak";
      }else if (settingVO_.getType().equalsIgnoreCase("relative to total amount")){
        subTitle = "relative to the total of all quantified peaks";
      }else if (settingVO_.getType().equalsIgnoreCase("amount end-volume")){
        subTitle = "amount in end volume";
      }else if (settingVO_.getType().equalsIgnoreCase("conc. end-volume")){
        subTitle = "concentration in end volume";
      }else if (settingVO_.getType().equalsIgnoreCase("weight end-volume")){
        subTitle = "weight in end volume";  
      }else if (settingVO_.getType().equalsIgnoreCase("amount sample-volume")){
        subTitle = "amount in sample volume";  
      }else if (settingVO_.getType().equalsIgnoreCase("conc. sample-volume")){
        subTitle = "concentration in probe volume";
      }else if (settingVO_.getType().equalsIgnoreCase("weight sample-volume")){
        subTitle = "weight in sample volume";  
      }else if (settingVO_.getType().equalsIgnoreCase("relative to sample weight")){
        subTitle = "relative to sample weight";
      }else if (settingVO_.getType().equalsIgnoreCase("relation to protein content")){
        subTitle = "relative to protein content";
      }else if (settingVO_.getType().equalsIgnoreCase("relation to neutral lipid content")){
        subTitle = "relative to indepently measured neutrel lipid content";        
      }else if (settingVO_.getType().equalsIgnoreCase("relation to measured neutral lipid")){
        subTitle = "relative to total class content measured by MS";          
      }
    }
    
    barChart_ = new BarChart(doubleSided_.isSelected(),logarithmicScale_.isSelected(),valuesToPaint_,sdsToPaint_,groupName_,subTitle,standardText,this.moleculeNames_,
        useSampleLookup_ ? sampleLookup_ : null, originalValueNames_,lowestValue,highestValue,/*amountOfGroupValues,*/yAxisText_, sdText, colorChooser_,colorType_);
    barChartPanel_.add(barChart_,BorderLayout.CENTER);
    barChartPanel_.invalidate();
    barChartPanel_.updateUI();
  }
  
  @SuppressWarnings( { "unchecked", "rawtypes" } )
  private Vector extractLowestHighestValue(String valueGetterMethod,String sdGetterMethod, String seGetterMethod, int maxIsotopes, ResultDisplaySettingsVO settingVO, boolean absolute,String preferredUnit){
    Vector results = new Vector();
    Double highestTotalValue = 0d;
    Double lowestTotalValue = Double.MAX_VALUE;
    Hashtable<String,Hashtable<String,Double>> valuesToPaint = new Hashtable<String,Hashtable<String,Double>>();
    Hashtable<String,Hashtable<String,Double>> sdsToPaint = new Hashtable<String,Hashtable<String,Double>>();
    Hashtable<String,Hashtable<String,Double>> originalFileValues = new Hashtable<String,Hashtable<String,Double>>();
    try{
      Class[] classes = new Class[2];      
      classes[0] = int.class;
      classes[1] = ResultDisplaySettingsVO.class;
      Method valueGetter = ResultCompVO.class.getMethod(valueGetterMethod, classes);
      Method sdGetter = ResultCompGroupVO.class.getMethod(sdGetterMethod, classes);
      Method seGetter = ResultCompGroupVO.class.getMethod(seGetterMethod, classes);
      Object[] getterParameters = new Object[2];
      getterParameters[1] = settingVO;
      for (String molName : analysisResultsHash_.keySet()){
        Hashtable<String,Double> valuesToPaintMol = new Hashtable<String,Double>();
        Hashtable<String,Double> sdsToPaintMol = new Hashtable<String,Double>();
        for (String name : analysisResultsHash_.get(molName).keySet()){
          ResultCompVO compVO = analysisResultsHash_.get(molName).get(name);
          int isotopes = extractCorrectIsotopesSetting(maxIsotopes,settingVO,valueGetterMethod,compVO);
          if (isotopes<maxIsotopes)
            maxIsotopes_.setSelectedItem(String.valueOf(isotopes));
          getterParameters[0] = isotopes;
          double standValue = (Double)valueGetter.invoke(compVO, getterParameters);
          if (absolute)
            standValue = StaticUtils.getAreaInCorrespondingUnit(standValue,preferredUnit);
          valuesToPaintMol.put(name, standValue);
          if (standValue>highestTotalValue)
            highestTotalValue = standValue;
          if (standValue>0&&standValue<lowestTotalValue)
            lowestTotalValue = standValue;
          if (type_ != TYPE_EXPERIMENT){
            if (!originalFileValues.containsKey(molName))
              originalFileValues.put(molName, new Hashtable<String,Double>());
          }else{
            if (!originalFileValues.containsKey(name))
              originalFileValues.put(name, new Hashtable<String,Double>());            
          }
          if (compVO instanceof ResultCompGroupVO){
            ResultCompGroupVO groupVO = (ResultCompGroupVO)compVO;
            double standValueSD = (Double)sdGetter.invoke(groupVO, getterParameters);
            if (absolute)
              standValueSD = StaticUtils.getAreaInCorrespondingUnit(standValueSD,preferredUnit);
            if (this.standardDeviation_!=null){
              if (this.standardDeviation_.isSelected()){
                standValueSD = standValueSD* Double.parseDouble(multStandardDeviation_.getText());
              }else{
                standValueSD = (Double)seGetter.invoke(groupVO, getterParameters);
                if (absolute)
                  standValueSD = StaticUtils.getAreaInCorrespondingUnit(standValueSD,preferredUnit);                
              }
            }
            sdsToPaintMol.put(name, standValueSD);
            if (Double.valueOf(standValueSD)!=Double.NaN && Double.valueOf(standValueSD)!=Double.NEGATIVE_INFINITY &&
                Double.valueOf(standValueSD)!=Double.POSITIVE_INFINITY){
              if ((standValue+standValueSD)>highestTotalValue)
                highestTotalValue = standValue+standValueSD;
              if ((standValue-standValueSD)>0&&(standValue-standValueSD)<lowestTotalValue)
                lowestTotalValue = standValue-standValueSD;
            }
            for (String expName : groupVO.getGroupingPartners().keySet()){
              if (type_ != TYPE_EXPERIMENT){
                if (originalFileValues.get(molName).containsKey(expName))
                  continue;
              }else{
                if (originalFileValues.get(name).containsKey(expName))
                  continue;                
              }
              ResultCompVO expCompVO = groupVO.getGroupingPartners().get(expName);
              double expStandValue = (Double)valueGetter.invoke(expCompVO, getterParameters);
              if (absolute)
                expStandValue = StaticUtils.getAreaInCorrespondingUnit(expStandValue,preferredUnit);
              if (type_ != TYPE_EXPERIMENT)
                originalFileValues.get(molName).put(expName, expStandValue);
              else
                originalFileValues.get(name).put(expName, expStandValue);
            }
          }else{
            if (type_ != TYPE_EXPERIMENT)
              originalFileValues.get(molName).put(name, standValue);
            else
              originalFileValues.get(name).put(name, standValue);
          }
        }
        valuesToPaint.put(molName, valuesToPaintMol);
        sdsToPaint.put(molName, sdsToPaintMol);
      }
    } catch(Exception ex){
      ex.printStackTrace();
    }
    results.add(lowestTotalValue);
    results.add(highestTotalValue);
    results.add(valuesToPaint);
    results.add(sdsToPaint);
    results.add(originalFileValues);
    return results;
  }
  
  private boolean verifySDFactInput(){
    try{
      Double.parseDouble(multStandardDeviation_.getText());
      return true;
    } catch (NumberFormatException nfx){
      new WarningMessage(new JFrame(), "Error", "SD value input must be in double format (xxx.xxx) and not "+multStandardDeviation_.getText());
      return false;
    }    
  }
  
  
  private class ReturnKeyListener implements java.awt.event.KeyListener{

    public void keyPressed(KeyEvent e)
    {
    }

    public void keyReleased(KeyEvent e)
    {
      if (e.getKeyCode() == KeyEvent.VK_ENTER){
        if (verifySDFactInput())
          createPaintingArea();
      }
    }

    public void keyTyped(KeyEvent e)
    {     
    }
    
  }
  
  private class DoubleVerifier extends InputVerifier{

    public boolean verify(JComponent input)
    {
      return verifySDFactInput();
    }
    
  }
  
  private class SDFocusListener implements FocusListener{

    public void focusGained(FocusEvent e)
    {
    }

    public void focusLost(FocusEvent e)
    {
      createPaintingArea();
      
    }
    
  }
  
  private class SelectionItemListener implements java.awt.event.ItemListener
  {
    private String m_ctrl;
    
    public SelectionItemListener(String ctrl){
      m_ctrl = ctrl;
    }

    public void itemStateChanged(ItemEvent e)  {
      if (m_ctrl.equalsIgnoreCase("ChangeChart")){
        if (e.getStateChange()==ItemEvent.SELECTED){
          repaintingBarChart_ = true;
          if (((String)quantType_.getSelectedItem()).equalsIgnoreCase("area absolute")){
            singleSided_.setSelected(true);
            singleSided_.setEnabled(false);
            doubleSided_.setEnabled(false);
            disableMagnitudeChooser();
            createPaintingArea();
          }
          if (((String)quantType_.getSelectedItem()).equalsIgnoreCase("percentual values")){
            singleSided_.setSelected(true);
            singleSided_.setEnabled(false);
            doubleSided_.setEnabled(false);
            disableMagnitudeChooser();
            createPaintingArea();
          }
          if (((String)quantType_.getSelectedItem()).equalsIgnoreCase("area relative to molecule")){
            doubleSided_.setSelected(true);
            singleSided_.setEnabled(true);
            doubleSided_.setEnabled(true);
            disableMagnitudeChooser();
            createPaintingArea();
          }
          if (((String)quantType_.getSelectedItem()).equalsIgnoreCase("area relative to standard")){
            doubleSided_.setSelected(true);
            singleSided_.setEnabled(true);
            doubleSided_.setEnabled(true);
            disableMagnitudeChooser();
            createPaintingArea();
          }
          if (((String)quantType_.getSelectedItem()).equalsIgnoreCase("absolute quantity")){
            singleSided_.setSelected(true);
            singleSided_.setEnabled(false);
            doubleSided_.setEnabled(false);
            enableMagnitudeChooser();
            createPaintingArea();
          }
          repaintingBarChart_ = false;
        }
      }
      if (m_ctrl.equalsIgnoreCase("ChangeSided")||m_ctrl.equalsIgnoreCase("ChangeLog")||m_ctrl.equalsIgnoreCase("ChangeMedian")){
        if (m_ctrl.equalsIgnoreCase("ChangeSided")){
          if (singleSided_.isSelected()){
            logarithmicScale_.setText("logarithmic");
            logarithmicScale_.setToolTipText(TooltipTexts.BARCHART_LOGARITHMIC);
            normalScale_.setText("linear");
            normalScale_.setToolTipText(TooltipTexts.BARCHART_LINEAR);
          }else{
            logarithmicScale_.setText("log10");
            logarithmicScale_.setToolTipText(TooltipTexts.BARCHART_LOG10);
            normalScale_.setText("log2");           
            normalScale_.setToolTipText(TooltipTexts.BARCHART_LOG2);
          }
        }
        if (!repaintingBarChart_){
          if (e.getStateChange()==ItemEvent.SELECTED){       
            createPaintingArea();
          }
        }
      }
      if (m_ctrl.equalsIgnoreCase("ChangeSD")){
        if (e.getStateChange()==ItemEvent.SELECTED){
          if (standardDeviation_.isSelected())
            multStandardDeviation_.setEnabled(true);
          else
            multStandardDeviation_.setEnabled(false);
          createPaintingArea();
        }
      }
      if (m_ctrl.equalsIgnoreCase("ChangeIsotope")){
        if (e.getStateChange()==ItemEvent.SELECTED){
          createPaintingArea();
        }
      }
      if (m_ctrl.equalsIgnoreCase("ChangeMagnitude")){
        if (!repaintingBarChart_){
          createPaintingArea();
        }
      }
    }
  }
  

  @SuppressWarnings("unchecked")
  public void actionPerformed(ActionEvent e)
  {
    if (e.getActionCommand().equalsIgnoreCase(ExportSettingsPanel.CHANGE_RT_SELECTION_STATUS)){
      if (exportRTDev_!=null){
        if (exportRT_!=null && exportRT_.isSelected()){
          exportRTDev_.setEnabled(true);
        }else{
          exportRTDev_.setSelected(false);
          exportRTDev_.setEnabled(false);
        }
      }
    }else if (e.getActionCommand().equalsIgnoreCase(ExportPanel.EXPORT_PNG)){
      exportFileChooser_.setFileFilter(new FileNameExtensionFilter("PNG (*.png)","png"));
      int returnVal = exportFileChooser_.showSaveDialog(new JFrame());
      if (returnVal == JFileChooser.APPROVE_OPTION) {
        File fileToStore = exportFileChooser_.getSelectedFile();
        @SuppressWarnings("rawtypes")
        Vector results = StaticUtils.checkFileStorage(fileToStore,"png",this);
        fileToStore = (File)results.get(0);
        if ((Boolean)results.get(1)){
          try {
            ImageIO.write(barChart_.getImage(), "PNG", fileToStore);
          }catch (IOException ex){ new WarningMessage(new JFrame(), "Error", ex.getMessage());}
        }  
      }
    }else if (e.getActionCommand().equalsIgnoreCase(ExportPanel.EXPORT_SVG)){
      exportFileChooser_.setFileFilter(new FileNameExtensionFilter("SVG (*.svg)","svg"));
      int returnVal = exportFileChooser_.showSaveDialog(new JFrame());
      if (returnVal == JFileChooser.APPROVE_OPTION) {
        File fileToStore = exportFileChooser_.getSelectedFile();
        @SuppressWarnings("rawtypes")
        Vector results = StaticUtils.checkFileStorage(fileToStore,"svg",this);
        fileToStore = (File)results.get(0);
        if ((Boolean)results.get(1)){
          try {
            DOMImplementation domImpl = GenericDOMImplementation.getDOMImplementation();

            // Create an instance of org.w3c.dom.Document.
            String svgNS = "http://www.w3.org/2000/svg";
            Document document = domImpl.createDocument(svgNS, "svg", null);

            // Create an instance of the SVG Generator.
            SVGGraphics2D svgGenerator = new SVGGraphics2D(document);
            barChart_.drawBarchart(svgGenerator,true);
            boolean useCSS = true; // we want to use CSS style attributes
            BufferedOutputStream stream = new BufferedOutputStream(new FileOutputStream(fileToStore));
            Writer out = new OutputStreamWriter(stream, "UTF-8");
            svgGenerator.stream(out, useCSS);
            out.close();
          }catch (IOException ex){ new WarningMessage(new JFrame(), "Error", ex.getMessage());}
        }  
      }      
    } else if (e.getActionCommand().equalsIgnoreCase(ExportPanel.EXPORT_EXCEL)||
        e.getActionCommand().equalsIgnoreCase(ExportPanel.EXPORT_TEXT)){
      int exportType =  ExportOptionsVO.EXPORT_NO_DEVIATION;
      String sdValue = null;
      if (this.standardDeviation_!=null  && this.exportDeviation_.isSelected()){
        if (standardDeviation_.isSelected()){
          exportType =  ExportOptionsVO.EXPORT_SD_DEVIATION;
          sdValue = multStandardDeviation_.getText();
        }else{
          exportType =  ExportOptionsVO.EXPORT_SD_ERROR;
        }
      }
      boolean exportRTDev = false;
      if (exportRTDev_!=null && exportRTDev_.isSelected())
        exportRTDev = true;
      boolean exportRT = false;
      if (exportRT_!=null && exportRT_.isSelected())
        exportRT = true;
      int roundValue = 6;
      if ((barChart_.getRoundValue()+2)>roundValue)
        roundValue = barChart_.getRoundValue()+2;
      //TODO: add radio buttons for the species type
      ExportOptionsVO exportVO = new ExportOptionsVO(exportType,sdValue,columnAnalyte_.isSelected(),exportRT,exportRTDev,roundValue,
          LipidomicsConstants.EXPORT_ANALYTE_TYPE_SPECIES);
      if (e.getActionCommand().equalsIgnoreCase(ExportPanel.EXPORT_EXCEL)){
        exportFileChooser_.setFileFilter(new FileNameExtensionFilter("Microsoft Office Excel Woorkbook (*.xlsx)","xlsx"));
        int returnVal = exportFileChooser_.showSaveDialog(new JFrame());
        if (returnVal == JFileChooser.APPROVE_OPTION) {
          File fileToStore = exportFileChooser_.getSelectedFile();
          @SuppressWarnings("rawtypes")
          Vector results = StaticUtils.checkFileStorage(fileToStore,"xlsx",this);
          fileToStore = (File)results.get(0);
          if ((Boolean)results.get(1)){
            try {
              int maxIsotope = Integer.parseInt((String)maxIsotopes_.getSelectedItem());
              @SuppressWarnings("rawtypes")
              Vector resultsVector = extractResultsVector();
              String headerTitle = StaticUtils.getAreaTypeString(settingVO_);
              ExcelAndTextExporter.exportToFile(type_!=TYPE_CLASSES,exportVO.getSpeciesType(), groupName_,new BufferedOutputStream(new FileOutputStream(fileToStore)),true,
                  maxIsotope, (Vector<String>)resultsVector.get(0), rtGrouped_, isGroupedView_, (Vector<String>)resultsVector.get(1),
                  (Hashtable<String,String>)resultsVector.get(2), sampleLookup_.getSampleResultFullPaths(), sampleLookup_.getSamplesOfGroups(),
                  (Hashtable<String,Hashtable<String,Vector<Double>>>)resultsVector.get(3), yAxisText_, headerTitle, exportVO, sampleLookup_.getComparativeResultsLookup(),
                  modifications_);

            }
            catch (NumberFormatException ex) {ex.printStackTrace();new WarningMessage(new JFrame(), "Error", ex.getMessage());}
            catch (FileNotFoundException ex) {ex.printStackTrace();new WarningMessage(new JFrame(), "Error", ex.getMessage());}
            catch (IOException ex) {ex.printStackTrace();new WarningMessage(new JFrame(), "Error", ex.getMessage());}
            catch (ExcelInputFileException | ExportException | SpectrummillParserException ex) {
              ex.printStackTrace();
              new WarningMessage(new JFrame(), "Error", ex.getMessage());
            }
          }  
        }
      } else if (e.getActionCommand().equalsIgnoreCase(ExportPanel.EXPORT_TEXT)){
        exportFileChooser_.setFileFilter(new FileNameExtensionFilter("Text (*.txt)","txt"));
        int returnVal = exportFileChooser_.showSaveDialog(new JFrame());
        if (returnVal == JFileChooser.APPROVE_OPTION) {
          File fileToStore = exportFileChooser_.getSelectedFile();
          @SuppressWarnings("rawtypes")
          Vector results = StaticUtils.checkFileStorage(fileToStore,"txt",this);
          fileToStore = (File)results.get(0);
          if ((Boolean)results.get(1)){
            try {
              int maxIsotope = Integer.parseInt((String)maxIsotopes_.getSelectedItem());
              @SuppressWarnings("rawtypes")
              Vector resultsVector = extractResultsVector();
              String headerTitle = StaticUtils.getAreaTypeString(settingVO_);
                if (((String)quantType_.getSelectedItem()).equalsIgnoreCase("area absolute")){
                  headerTitle = "Area";
                } else if (((String)quantType_.getSelectedItem()).equalsIgnoreCase("area relative to molecule")){
                  headerTitle = "Relative to Median";
                } else if (((String)quantType_.getSelectedItem()).equalsIgnoreCase("area relative to standard")){
                  headerTitle = "Relative to Standard";
                }
                ExcelAndTextExporter.exportToFile(type_!=TYPE_CLASSES,exportVO.getSpeciesType(), groupName_,new BufferedOutputStream(new FileOutputStream(fileToStore)),false,
                    maxIsotope, (Vector<String>)resultsVector.get(0), rtGrouped_, isGroupedView_, (Vector<String>)resultsVector.get(1),
                    (Hashtable<String,String>)resultsVector.get(2), sampleLookup_.getSampleResultFullPaths(), sampleLookup_.getSamplesOfGroups(),
                    (Hashtable<String,Hashtable<String,Vector<Double>>>)resultsVector.get(3), yAxisText_, headerTitle, exportVO, sampleLookup_.getComparativeResultsLookup(),
                    modifications_);

            }
            catch (NumberFormatException ex) {new WarningMessage(new JFrame(), "Error", ex.getMessage());}
            catch (FileNotFoundException ex) {new WarningMessage(new JFrame(), "Error", ex.getMessage());}
            catch (IOException ex) {new WarningMessage(new JFrame(), "Error", ex.getMessage());}
            catch (ExcelInputFileException | ExportException | SpectrummillParserException ex) {
              new WarningMessage(new JFrame(), "Error", ex.getMessage());
            }
          }  
        }      
      }
    } 
  }
  
  
  @SuppressWarnings({ "unchecked", "rawtypes" })
  private Vector extractResultsVector(){
    Vector results = new  Vector();
    Vector<String> molNames = new Vector<String>();
    Vector<String> expNames = new Vector<String>();
    Hashtable<String,String> expHash = new Hashtable<String,String>();
    Hashtable<String,Hashtable<String,Vector<Double>>> resultsHash = new Hashtable<String,Hashtable<String,Vector<Double>>>();
    if (type_ == TYPE_MOLECULE || type_ == TYPE_CLASSES){
      molNames.addAll(moleculeNames_);
      expNames.addAll(originalValueNames_);
      for (String expName : originalValueNames_)
        expHash.put(expName, sampleLookup_.getDisplayName(expName));
      for (String molName :  this.moleculeNames_){
        Hashtable<String,Vector<Double>> resultHash = new Hashtable<String,Vector<Double>>();
        for (String expName : sampleLookup_.getComparativeResultsLookup().getExpNamesInSequence()){/////for (String expName : originalValueNames_){// old: BarChart.extractDisplayNames(sampleLookup_,originalValueNames_
          Vector<Double> resultsVect = new Vector<Double>();
          resultsVect.add(singleResultValues_.get(molName).get(expName));
          String displayName = expName;
          resultHash.put(displayName, resultsVect);
        }  
        resultsHash.put(molName, resultHash);
      }
    } else if (type_ == TYPE_EXPERIMENT){
      molNames.addAll(originalValueNames_);
      String rightExp = null;
      if (isGroupedView_)
        rightExp = moleculeNames_.get(0);
      else{
        for (String expName : sampleLookup_.getComparativeResultsLookup().getExpNamesInSequence()){
          if (!sampleLookup_.getDisplayName(expName).equalsIgnoreCase(moleculeNames_.get(0)))
            continue;
          rightExp = expName;
        }
      }
      expNames.add(rightExp);
      expHash.put(rightExp, moleculeNames_.get(0));
      for (String molName : valuesToPaint_.values().iterator().next().keySet()){
        Hashtable<String,Vector<Double>> resultHash = new Hashtable<String,Vector<Double>>();
        if (isGroupedView_){
          for (String expName : sampleLookup_.getSamplesOfGroups().get(rightExp)){
            Vector<Double> resultsVect = new Vector<Double>();
            resultsVect.add(singleResultValues_.get(molName).get(expName));
            resultHash.put(expName, resultsVect);
          }
        }else{
          Vector<Double> resultsVect = new Vector<Double>();
          resultsVect.add(valuesToPaint_.values().iterator().next().get(molName));
          resultHash.put(rightExp, resultsVect);
        }
        resultsHash.put(molName, resultHash);
      }  
      
    }
    results.add(molNames);
    results.add(expNames);
    
    results.add(expHash);
    results.add(resultsHash);
    return results;
  }
  
  private int extractCorrectIsotopesSetting(int maxIsotopes,ResultDisplaySettingsVO settingVO, String valueGetterMethod, ResultCompVO compVO){
    int isotopes = maxIsotopes;
    if ((!settingVO.getType().equalsIgnoreCase("relative to measured class amount") && 
        !settingVO.getType().equalsIgnoreCase("relative to total amount")) || valueGetterMethod.equalsIgnoreCase("getRelativeValue")){
      isotopes = compVO.getAvailableIsotopeNr(maxIsotopes);
    }  
    return isotopes;
  }
  
  private void calculatePercentualValues(int maxIsotope) throws CalculationNotPossibleException{
    ResultDisplaySettingsVO settingVO = new ResultDisplaySettingsVO(settingVO_);
    settingVO.setDivisorMagnitude("");
    boolean isGrouped = false;
    Vector<String> separateEntities_ = new Vector<String>();
    Vector<String> valuesToSum_ = new Vector<String>();
    if (bigConstructor_){
      valuesToSum_ = moleculeNames_;
      separateEntities_ = new Vector<String>(originalValueNames_);//BarChart.extractDisplayNames(sampleLookup_,originalValueNames_);
    }else{
      valuesToSum_= new Vector<String>(originalValueNames_);//BarChart.extractDisplayNames(sampleLookup_,originalValueNames_);
      separateEntities_ = moleculeNames_;
    }
    for (String valueName : separateEntities_){
      double totalArea = 0;
      double totalStd = 0;
      Hashtable<String,Double> singleExpAreas = new Hashtable<String,Double>();
      for (String molName : valuesToSum_){
        ResultCompVO compVO = null;
        if (bigConstructor_)
          compVO = analysisResultsHash_.get(molName).get(valueName);
        else
          compVO = analysisResultsHash_.get(valueName).get(molName);
        int isotopes = this.extractCorrectIsotopesSetting(maxIsotope, settingVO, "", compVO);
        double area = compVO.getArea(isotopes, settingVO);
        if (area>0 && !Double.isNaN(area) && !Double.isInfinite(area))
          totalArea += area;
        if (compVO instanceof ResultCompGroupVO){
          isGrouped = true;
          ResultCompGroupVO groupVO = (ResultCompGroupVO)compVO;
          for (String origExpName : groupVO.getGroupingPartners().keySet()){
            double singleExpArea = 0;
            if (singleExpAreas.containsKey(origExpName)) singleExpArea = singleExpAreas.get(origExpName);
            singleExpArea += groupVO.getGroupingPartners().get(origExpName).getArea(compVO.getAvailableIsotopeNr(isotopes),settingVO);
            singleExpAreas.put(origExpName, singleExpArea);
          }
        }
      }

      Vector<Double> sumAreas = new Vector<Double>();
      Hashtable<String,Vector<Double>> singleAreas = new Hashtable<String,Vector<Double>>();
      for (int i=0; i!=maximumIsotopes_;i++){
        sumAreas.add(totalArea);
        for (String origExp :singleExpAreas.keySet()){
          Vector<Double> singleArea = new Vector<Double>();
          if (singleAreas.containsKey(origExp)) singleArea =  singleAreas.get(origExp);
          //if single analytes are compared, no normalisation is required; just the analytes have to be summed as they are
          if (!bigConstructor_ && type_ == TYPE_MOLECULE){
            singleArea.add(1d);
          }else
            singleArea.add(singleExpAreas.get(origExp));
          singleAreas.put(origExp,singleArea);
        }
      }

      Vector<Double> sumMeans = new Vector<Double>();
      Vector<Double> sumSds = new Vector<Double>();
      if (isGrouped){
        double totalMean = 0;
        Vector<Double> stdevs = new Vector<Double>();
        for (String molName : valuesToSum_){
          ResultCompGroupVO compVO = null;
          if (bigConstructor_)
            compVO = (ResultCompGroupVO)analysisResultsHash_.get(molName).get(valueName);
          else
            compVO = (ResultCompGroupVO)analysisResultsHash_.get(valueName).get(molName);
          Vector<Double> valuesForMean = new Vector<Double>();
          int isotopes = this.extractCorrectIsotopesSetting(maxIsotope, settingVO, "", compVO);
          for (String origExpName : compVO.getGroupingPartners().keySet()){
            ResultCompVO groupVO = compVO.getGroupingPartners().get(origExpName);
            groupVO.setSumValueForPercentage(singleAreas.get(origExpName ));
            ResultDisplaySettingsVO newVO = new ResultDisplaySettingsVO(settingVO);
            newVO.setPercent(true);
            newVO.setDivisorMagnitude("");
            double ratio = groupVO.getRatioToPercentualValue(groupVO.getAvailableIsotopeNr(isotopes),newVO);
            if (ratio>0 && !Double.isNaN(ratio) && !Double.isInfinite(ratio))
              valuesForMean.add(ratio);
          }
          double mean = Calculator.mean(valuesForMean);
          double stdev = Calculator.stddeviation(valuesForMean);
          if (mean>0 && !Double.isNaN(mean) && !Double.isInfinite(mean)){
            totalMean+=mean;
            stdevs.add(stdev);
          }
        }     
        totalStd = Calculator.calculateSumStdevErrorPropagated(stdevs);
        for (int i=0; i!=maximumIsotopes_;i++){
          sumMeans.add(totalMean);
          sumSds.add(totalStd);
        }
      }
      for (String molName : valuesToSum_){
        ResultCompVO compVO = null;
        if (bigConstructor_)
          compVO = analysisResultsHash_.get(molName).get(valueName);
        else
          compVO = analysisResultsHash_.get(valueName).get(molName);
        compVO.setSumValueForPercentage(sumAreas);
        if (isGrouped){
          ((ResultCompGroupVO)compVO).setSumValueForPercentage(sumAreas);
          ((ResultCompGroupVO)compVO).setSumPercentualMeans(sumMeans);
          ((ResultCompGroupVO)compVO).setSumPercentualSds(sumSds);
          Hashtable<String,ResultCompVO> hashes = ((ResultCompGroupVO)compVO).getGroupingPartners();
          for (String origExp : hashes.keySet()){
            ResultCompVO groupVO = hashes.get(origExp);
            
            groupVO.setSumValueForPercentage(singleAreas.get(origExp));
          }

        }  
      }
    }
  }
  
  private String getStandardName(Hashtable<String,Integer> hash, int isMethod){
    String standName = "";
    for (String oneStandName : hash.keySet()){
      if (isMethod == hash.get(oneStandName))
        return oneStandName;
    }
    return standName;
  }
  
  private void addEmptyItemToChooser(){
    magnitudeChooser_.addItem("    ");
  }
  
  private void enableMagnitudeChooser(){
    magnitudeChooser_.removeAllItems();
    if (((String)quantType_.getSelectedItem()).equalsIgnoreCase("absolute quantity")){
      if (preferredUnit_.equalsIgnoreCase("")||preferredUnit_.equalsIgnoreCase("m")||
          preferredUnit_.equalsIgnoreCase("\u03BC")||preferredUnit_.equalsIgnoreCase("n")||
          preferredUnit_.equalsIgnoreCase("p")||preferredUnit_.equalsIgnoreCase("f")||
          preferredUnit_.equalsIgnoreCase("a")){
        magnitudeChooser_.addItem("");
        magnitudeChooser_.addItem("m");
        magnitudeChooser_.addItem("\u03BC");
        magnitudeChooser_.addItem("n");
        magnitudeChooser_.addItem("p");
        magnitudeChooser_.addItem("f");
        magnitudeChooser_.addItem("a");
      }else if (preferredUnit_.equalsIgnoreCase("%")|| preferredUnit_.equalsIgnoreCase("\u2030")){
        magnitudeChooser_.addItem("%");
        magnitudeChooser_.addItem("\u2030");
      }
      magnitudeChooser_.setSelectedItem(preferredUnit_);
    }
    magnitudeChooser_.setEnabled(true);
  }
  
  private void disableMagnitudeChooser(){
    magnitudeChooser_.removeAllItems();
    addEmptyItemToChooser();
    magnitudeChooser_.setEnabled(false);
  }
  
}
