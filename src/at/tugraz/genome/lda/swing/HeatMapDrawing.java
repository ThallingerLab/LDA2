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
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.MenuItem;
import java.awt.Point;
import java.awt.PopupMenu;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.MouseEvent;
import java.awt.image.BufferedImage;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
import java.util.List;
import java.util.Timer;
import java.util.TimerTask;
import java.util.Vector;

import javax.imageio.ImageIO;
import javax.swing.Icon;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.SwingUtilities;
import javax.swing.event.MouseInputAdapter;
import javax.swing.filechooser.FileNameExtensionFilter;

import org.apache.batik.dom.GenericDOMImplementation;
import org.apache.batik.svggen.SVGGraphics2D;
import org.w3c.dom.DOMImplementation;
import org.w3c.dom.Document;

import at.tugraz.genome.lda.ChromExportThread;
import at.tugraz.genome.lda.LipidDataAnalyzer;
import at.tugraz.genome.lda.TooltipTexts;
import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.analysis.HeatMapClickListener;
import at.tugraz.genome.lda.analysis.LipidomicsHeatMap;
import at.tugraz.genome.lda.analysis.exception.CalculationNotPossibleException;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.AutoAnalyteAddVO;
import at.tugraz.genome.lda.vos.ExportOptionsVO;
import at.tugraz.genome.lda.vos.ResultCompGroupVO;
import at.tugraz.genome.lda.vos.ResultCompVO;
import at.tugraz.genome.lda.vos.ResultDisplaySettingsVO;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class HeatMapDrawing extends JPanel implements ActionListener
{

  /**
   * 
   */
  private static final long serialVersionUID = 9212631108743680488L;
  protected static final String CHANGE_IS_STATUS = "changeISShowStatus";
  protected static final String CHANGE_ES_STATUS = "changeESShowStatus";
  protected static final String CHANGE_DOUBLE_STATUS = "changeDoublePeaksStatus";
  
  private LipidomicsHeatMap heatmap_;
  private BufferedImage renderedImage_;
  private int imagePositionY_ = 0;
  private int imagePositionX_ = 0;
  private Rectangle rectToDraw_ = null;
  private JLabel label_;
  private JPanel heatMapPanel_;
  private JLabel statusText_;
  private JPanel exportProgressPanel_;
  protected ResultCompVO[][] compVOs_;
  protected HeatMapClickListener heatMapListener_;
  protected String groupName_;
  protected JCheckBox showInternalStandards_;
  protected JCheckBox showExternalStandards_;
  protected JCheckBox markDoublePeaks_;
  protected JComboBox<String> maxIsotopes_;
  protected Hashtable<String,Hashtable<String,ResultCompVO>> resultsOfOneGroup_;
  protected Vector<String> experimentNames_;
  protected Vector<String> moleculeNames_;
  protected Hashtable<String,Integer> isLookup_;
  protected Hashtable<String,Integer> esLookup_;
  protected ResultDisplaySettingsVO settingsVO_;
  protected ResultDisplaySettings displaySettings_;
  protected ResultSelectionSettings selectionSettings_;
  protected ExportSettingsPanel exportSettings_;
  protected ResultSelectionSettings combinedChartSettings_;
  protected boolean isGrouped_ = false;
  protected Hashtable<String,String> preferredUnit_ = null;
  private JFileChooser exportFileChooser_;
  private String molGroupName_;
  protected JLabel exportLabel_;
  protected JProgressBar exportProgress_;
  protected JLabel spinnerLabel_;
  protected JButton cancelExport_;
  
  
  private Hashtable<String,Hashtable<String,Color>> attentionProbes_;
  private Hashtable<String,String> selectedMolecules_;
  
  private boolean combinedDialogOpen_ = false;
  private boolean parentAction_;
  
  private PopupMenu applySettingsPopup_;
  private PopupMenu removeAnalytePopup_;
  
  private MenuItem selectItem_;
  private MenuItem deselectItem_;
  
  private MenuItem applyToAllDoubles_;
//  private String lastClickedExp_;
//  private String lastClickedMol_;
//  private ResultCompVO lastClickedResultVO_;
  private int[] lastClickedCellPos_;
  private String lastClickedAnalyte_;
  private Vector<String> modifications_;
  private ChromExportDialog chromExport_;
  private Hashtable<String,String> fromShortToExpName_;
  private Double rtTolerance_;
  
  private ChromExportThread chromExportThread_;
  private Timer timer_;

  public HeatMapDrawing(String molGroupName,Hashtable<String,Hashtable<String,ResultCompVO>> resultsOfOneGroup, Vector<String> experimentNames, 
      Vector<String> moleculeNames, Hashtable<String,Integer> isLookup, Hashtable<String,Integer> esLookup, int maxIsotopesOfGroup,
      Vector<String> modifications, JLabel statusText, HeatMapClickListener listener, String groupName, boolean isGrouped,
      boolean isAvailability, boolean esAvailability, boolean absoluteSettings, boolean hasProtein,
      boolean hasNeutralLipid, ResultDisplaySettings displaySettings, ResultSelectionSettings selectionSettings,
      ResultSelectionSettings combinedChartSettings, Double rtTolerance){
    selectedMolecules_ = new Hashtable<String,String>();
    parentAction_ = true;
    this.rtTolerance_ = rtTolerance;
    molGroupName_ = molGroupName;
    combinedDialogOpen_ = false;
    this.isGrouped_ = isGrouped;
    modifications_ = modifications;
    this.setLayout(new BorderLayout());
    JPanel bottomPanel = new JPanel();
    bottomPanel.setLayout(new GridBagLayout());
    this.add(bottomPanel,BorderLayout.SOUTH);
    showInternalStandards_  = new JCheckBox("show intern. stand.");
    showInternalStandards_.setSelected(true);
    showInternalStandards_.setActionCommand(CHANGE_IS_STATUS);
    showInternalStandards_.addActionListener(this);
    showInternalStandards_.setToolTipText(TooltipTexts.HEATMAP_SHOW_INT);
    bottomPanel.add(showInternalStandards_, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 2, 2, 2), 0, 0));
    
    maxIsotopes_ = new JComboBox<String>();
    for (int i=0;i!=maxIsotopesOfGroup;i++){
      maxIsotopes_.addItem(String.valueOf(i));
    }
    maxIsotopes_.setSelectedItem(String.valueOf(maxIsotopesOfGroup-1));
    SelectionItemListener isotopeListener = new SelectionItemListener("ChangeIsotope");
    maxIsotopes_.addItemListener(isotopeListener);
    maxIsotopes_.setToolTipText(TooltipTexts.HEATMAP_ISOTOPES);
    bottomPanel.add(maxIsotopes_, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 2, 2, 2), 0, 0));
    bottomPanel.add(new JLabel("isotopes"), new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 2, 2, 2), 0, 0));

    showExternalStandards_  = new JCheckBox("show extern. stand.");
    showExternalStandards_.setSelected(true);
    showExternalStandards_.setActionCommand(CHANGE_ES_STATUS);
    showExternalStandards_.addActionListener(this);
    showExternalStandards_.setToolTipText(TooltipTexts.HEATMAP_SHOW_EXT);
    bottomPanel.add(showExternalStandards_, new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 2, 2, 2), 0, 0));

    if (!isGrouped_){
      markDoublePeaks_  = new JCheckBox("double peaks/missed mods");
      markDoublePeaks_.setSelected(true);
      markDoublePeaks_.setActionCommand(CHANGE_DOUBLE_STATUS);
      markDoublePeaks_.addActionListener(this);
      markDoublePeaks_.setToolTipText(TooltipTexts.HEATMAP_DOUBLE_PEAKS);
      bottomPanel.add(markDoublePeaks_, new GridBagConstraints(0, 2, 2, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 2, 2, 2), 0, 0));
    }
    
    JPanel buttonsPanel = new JPanel();
    JButton settingsButton = new JButton("Settings");
    settingsButton.addActionListener(this);
    settingsButton.setMargin(new Insets(1,5,1,5));
    settingsButton.setActionCommand("openSettingsDialog");
    settingsButton.setToolTipText(TooltipTexts.HEATMAP_SETTINGS);
    buttonsPanel.add(settingsButton);

    JButton selectionButton = new JButton("Select molecules");
    selectionButton.addActionListener(this);
    selectionButton.setMargin(new Insets(1,5,1,5));
    selectionButton.setActionCommand("openSelectionDialog");
    selectionButton.setToolTipText(TooltipTexts.HEATMAP_SELECTED);
    buttonsPanel.add(selectionButton);
    
    bottomPanel.add(buttonsPanel, new GridBagConstraints(0, 3, 3, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(2, 2, 8, 2), 0, 0));
    
    JPanel buttonsPanel2 = new JPanel();
    JButton combinedButton = new JButton("Combined chart");
    combinedButton.addActionListener(this);
    combinedButton.setMargin(new Insets(1,5,1,5));
    combinedButton.setActionCommand("openCombinedDialog");
    combinedButton.setToolTipText(TooltipTexts.HEATMAP_COMBINED);
    buttonsPanel2.add(combinedButton);   
//    if (isGrouped){
    JButton exportButton = new JButton("Export options");
    exportButton.addActionListener(this);
    exportButton.setMargin(new Insets(1,5,1,5));
    exportButton.setActionCommand("exportSelectionDialog");
    exportButton.setToolTipText(TooltipTexts.HEATMAP_EXPORT_OPTIONS);
    buttonsPanel2.add(exportButton);
      
//    }
    bottomPanel.add(buttonsPanel2, new GridBagConstraints(0, 4, 3, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(2, 2, 8, 2), 0, 0));

    applySettingsPopup_ = new PopupMenu("Apply Peak-RT to double peaks");
    applyToAllDoubles_ = new MenuItem("Choose just one peak for doubles");
    applyToAllDoubles_.addActionListener(this);
    applySettingsPopup_.add(applyToAllDoubles_);
    MenuItem quantOthers = new MenuItem("Quant. anal. at not found");
    quantOthers.addActionListener(this);
    applySettingsPopup_.add(quantOthers);
    quantOthers = new MenuItem("Take exact peak for others");
    quantOthers.addActionListener(this);
    applySettingsPopup_.add(quantOthers);
    
    
    this .add(applySettingsPopup_);
    
    removeAnalytePopup_ = new PopupMenu("Remove analyte");
    MenuItem removeItem = new MenuItem("Remove analyte in all probes");
    removeItem.addActionListener(this);
    removeAnalytePopup_.add(removeItem);
    selectItem_ = new MenuItem("Select analyte");
    selectItem_.addActionListener(this);
    removeAnalytePopup_.add(selectItem_);
    deselectItem_ = new MenuItem("Deselect analyte");
    deselectItem_.addActionListener(this);
    removeAnalytePopup_.add(deselectItem_);    
    this .add(removeAnalytePopup_);
    
    this.resultsOfOneGroup_ = resultsOfOneGroup;
    this.experimentNames_ = new Vector<String>(experimentNames);
    this.moleculeNames_ = moleculeNames;
    this.isLookup_ = isLookup;
    this.esLookup_ = esLookup;
    this.groupName_ = groupName;
    statusText_ = statusText;
    this.heatMapListener_ = listener;
    displaySettings_ = displaySettings;
//    displaySettings_ = new ResultDisplaySettings(isAvailability,esAvailability,isLookup_,esLookup_,absoluteSettings, hasProtein,
//        hasNeutralLipid,this);
    settingsVO_ = displaySettings_.getSettingsVO();
    selectionSettings_ = selectionSettings;
    combinedChartSettings_ = combinedChartSettings;
//    selectionSettings_ = new ResultSelectionSettings(moleculeNames_,true,this);
//    combinedChartSettings_ = new ResultSelectionSettings(moleculeNames_,false,this);
    exportSettings_ = new ExportSettingsPanel(isGrouped_,this);
    
    exportFileChooser_ = new JFileChooser();
    exportFileChooser_.setPreferredSize(new Dimension(600,500));
    
    if (!isGrouped_){     
      chromExport_ = new ChromExportDialog("Chrom export",getDisplayNames(),moleculeNames_,modifications_,this);
////      chromExport_.addActionListener(this);
    }
    chromExportThread_ = null;
    this.initTimer();
    this.generateHeatMap();
  }
  
  @SuppressWarnings("unchecked")
  public void generateHeatMap(){
    preferredUnit_ = new Hashtable<String,String>();
    if (this.heatMapPanel_!=null)
      this.remove(heatMapPanel_);
    
    Vector<String> molNames = this.getSelectedMoleculeNames();
    
//    if (!this.showInternalStandards_.isSelected()){
//      Vector<String> molWithoutIs = new Vector<String>();
//      for (String name : molNames){
//        if (!isLookup_.containsKey(name))
//          molWithoutIs.add(name);
//      }
//      molNames = new Vector<String>(molWithoutIs);
//    }
    try{
//      if (settingsVO_.getType().equalsIgnoreCase("amount probe-volume")){
//        for (String expName : experimentNames_){
//          int maxIsotope = Integer.parseInt((String)maxIsotopes_.getSelectedItem());
//          double totalLipid = resultsOfOneGroup_.get(resultsOfOneGroup_.keySet().iterator().next()).get(expName).getTotalGroupMass(maxIsotope, settingsVO_.getISStandMethod(), settingsVO_.getESStandMethod());
//          System.out.println(expName+" ; "+(totalLipid*1000d)+"mg");
//        }
//      }
      int maxIsotopes = 0;
      try{
        maxIsotopes = Integer.parseInt((String)maxIsotopes_.getSelectedItem());
      } catch (NumberFormatException nfx){}  
      @SuppressWarnings("rawtypes")
      Vector results = LipidomicsHeatMap.getDataObjectForConstructor(resultsOfOneGroup_, experimentNames_, molNames,maxIsotopes, settingsVO_);
      attentionProbes_ = (Hashtable<String,Hashtable<String,Color>>)results.get(3);
      Hashtable<String,Hashtable<String,Color>> attentionProbesToPaint = new Hashtable<String,Hashtable<String,Color>>(attentionProbes_);
      if (markDoublePeaks_!=null && !markDoublePeaks_.isSelected())
        attentionProbesToPaint = new Hashtable<String,Hashtable<String,Color>>();
      this.heatmap_ = new LipidomicsHeatMap((float[][])results.get(0),experimentNames_,heatMapListener_, molNames, attentionProbesToPaint);
      this.compVOs_ = (ResultCompVO[][])results.get(1);
      preferredUnit_ = (Hashtable<String,String>)results.get(2);
      renderedImage_ = this.heatmap_.createImage();
    
      heatMapPanel_ = new JPanel();
      heatMapPanel_.setLayout(new GridBagLayout());
      JPanel justHeatMapPanel = new JPanel();
      label_ = new JLabel(new ImageIcon(renderedImage_)); 
      MyListener myListener = new MyListener();
      label_.addMouseListener(myListener);
      label_.addMouseMotionListener(myListener);
      justHeatMapPanel.add(label_);
      heatMapPanel_.add(justHeatMapPanel,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
      
      ExportPanel exportPanel = new ExportPanel(Color.BLACK,Color.WHITE,this,!isGrouped_,true, false);
      heatMapPanel_.add(exportPanel,new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
          ,GridBagConstraints.NORTHWEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
      exportProgressPanel_ = new JPanel();
      exportLabel_ = new JLabel("");
      exportLabel_.setToolTipText(TooltipTexts.EXPORT_STATUS_TEXT);
      exportProgressPanel_.add(exportLabel_);
      exportProgress_ = new JProgressBar();
      exportProgress_.setMaximum(100);
      exportProgress_.setToolTipText(TooltipTexts.EXPORT_PROGRESS);
      exportProgressPanel_.add(exportProgress_);
      Icon icon = new ImageIcon(LipidDataAnalyzer.class.getResource("/images/spinner.gif"));
      spinnerLabel_ = new JLabel(icon);
      exportProgressPanel_.add(spinnerLabel_);
      cancelExport_ = new JButton("Cancel");
      cancelExport_.addActionListener(this);
      cancelExport_.setActionCommand("stopChromExport");
      cancelExport_.setToolTipText(TooltipTexts.EXPORT_STOP);
      exportProgressPanel_.add(cancelExport_);
      
      //exportProgressPanel_.add(new JLabel("test"));
      exportProgressPanel_.setVisible(false);
      heatMapPanel_.add(exportProgressPanel_,new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
          ,GridBagConstraints.NORTHWEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
      this.add(heatMapPanel_,BorderLayout.NORTH);
    } catch (CalculationNotPossibleException cex){
      new WarningMessage(new JFrame(),"ERROR",cex.getMessage());
    }
    if (!isGrouped_) chromExport_.refreshNames(getDisplayNames());
  }
  
  
  public Dimension getTotalSize() {
    return new Dimension(renderedImage_.getWidth()+10,renderedImage_.getHeight());
  }
  
  public void actionPerformed(ActionEvent event){
    actionPerformed(event.getActionCommand());
  }
  
  private void actionPerformed(String actionCommand)
  {
    if (actionCommand.equalsIgnoreCase(CHANGE_IS_STATUS)||
        actionCommand.equalsIgnoreCase(CHANGE_ES_STATUS)||
        actionCommand.equalsIgnoreCase(CHANGE_DOUBLE_STATUS)||
        actionCommand.equalsIgnoreCase("AcceptDisplaySettings")){
      settingsVO_ = displaySettings_.getSettingsVO();
      selectionSettings_.setVisible(false);
      exportSettings_.setVisible(false);
      boolean update = true;    
      if (combinedDialogOpen_){
        update = false;
        if (combinedChartSettings_.getSelected().size()==0){
          new WarningMessage(new JFrame(), "Error", "You have to select at least one molecule!");
        }else{
          combinedChartSettings_.setVisible(false);
          combinedDialogOpen_ = false;
          int maxIsotopes = Integer.parseInt((String)maxIsotopes_.getSelectedItem());
          Hashtable<String,String> selectedHash = new Hashtable<String,String>();
          Hashtable<String,String> preferredUnits = new Hashtable<String,String>();
          for (String name : combinedChartSettings_.getSelected()) selectedHash.put(name, name);
          for (String name : preferredUnit_.keySet()){
            if (selectedHash.containsKey(name))
              preferredUnits.put(name, preferredUnit_.get(name));
          }
          String preferredUnit = extractPreferredUnitForExp(preferredUnits);
          if (isGrouped_)
            heatMapListener_.combinedAnalyteGroupSelected(combinedChartSettings_.getSelected(), groupName_, maxIsotopes, settingsVO_, preferredUnit, StaticUtils.getCorrespondingUnit(settingsVO_, preferredUnit,true));
          else
            heatMapListener_.combinedAnalyteSelected(combinedChartSettings_.getSelected(), groupName_, maxIsotopes, settingsVO_, preferredUnit, StaticUtils.getCorrespondingUnit(settingsVO_, preferredUnit,true));
        }
      }
      if (actionCommand.equalsIgnoreCase(CHANGE_IS_STATUS) && !isGrouped_){
        chromExport_.checkMolecules(getISs(),showInternalStandards_.isSelected());
      }else if (actionCommand.equalsIgnoreCase(CHANGE_ES_STATUS) && !isGrouped_){
        chromExport_.checkMolecules(getESs(),showExternalStandards_.isSelected());
      } else if (actionCommand.equalsIgnoreCase("AcceptDisplaySettings") && !isGrouped_ ){
        chromExport_.checkMolecules(getUnselectedMoleculeNames(),false);
      }
      if (actionCommand.equalsIgnoreCase(CHANGE_IS_STATUS) && parentAction_)
        heatMapListener_.changeISStatus(molGroupName_, isGrouped_, showInternalStandards_.isSelected());
      if (actionCommand.equalsIgnoreCase(CHANGE_ES_STATUS) && parentAction_)
        heatMapListener_.changeESStatus(molGroupName_, isGrouped_, showExternalStandards_.isSelected());
//      if (actionCommand.equalsIgnoreCase(CHANGE_DOUBLE_STATUS) && parentAction_)
//        heatMapListener_.changeDoublePeakStatus(molGroupName_, markDoublePeaks_.isSelected());
      if (update){
        this.generateHeatMap();
        this.invalidate();
        this.updateUI();
      }
    } else if (actionCommand.equalsIgnoreCase("openSettingsDialog")){
      displaySettings_.setVisible(true);
    } else if (actionCommand.equalsIgnoreCase("openSelectionDialog")){
      selectionSettings_.setVisible(true);
    } else if (actionCommand.equalsIgnoreCase("exportSelectionDialog")){
      exportSettings_.setVisible(true);
    } else if (actionCommand.equalsIgnoreCase("openCombinedDialog")){
      combinedChartSettings_.setVisible(true);
      combinedDialogOpen_ = true;
    } else if (actionCommand.equalsIgnoreCase(ExportPanel.EXPORT_PNG)){
      exportFileChooser_.setFileFilter(new FileNameExtensionFilter("PNG (*.png)","png"));
      int returnVal = exportFileChooser_.showSaveDialog(new JFrame());
      if (returnVal == JFileChooser.APPROVE_OPTION) {
        File fileToStore = exportFileChooser_.getSelectedFile();
        @SuppressWarnings("rawtypes")
        Vector results = this.checkFileStorage(fileToStore,"png");
        fileToStore = (File)results.get(0);
        if ((Boolean)results.get(1)){
          try {
            ImageIO.write(renderedImage_, "PNG", fileToStore);
          }catch (IOException e){ new WarningMessage(new JFrame(), "Error", e.getMessage());}
        }  
      }
    } else if (actionCommand.equalsIgnoreCase(ExportPanel.EXPORT_SVG)){
      exportFileChooser_.setFileFilter(new FileNameExtensionFilter("SVG (*.svg)","svg"));
      int returnVal = exportFileChooser_.showSaveDialog(new JFrame());
      if (returnVal == JFileChooser.APPROVE_OPTION) {
        File fileToStore = exportFileChooser_.getSelectedFile();
        @SuppressWarnings("rawtypes")
        Vector results = this.checkFileStorage(fileToStore,"svg");
        fileToStore = (File)results.get(0);
        if ((Boolean)results.get(1)){
          try {
            DOMImplementation domImpl = GenericDOMImplementation.getDOMImplementation();

            String svgNS = "http://www.w3.org/2000/svg";
            Document document = domImpl.createDocument(svgNS, "svg", null);

            SVGGraphics2D svgGenerator = new SVGGraphics2D(document);
            heatmap_.createImage(svgGenerator);
            boolean useCSS = true; // we want to use CSS style attributes
            BufferedOutputStream stream = new BufferedOutputStream(new FileOutputStream(fileToStore));
            Writer out = new OutputStreamWriter(stream, "UTF-8");
            svgGenerator.stream(out, useCSS);

          }catch (IOException e){ new WarningMessage(new JFrame(), "Error", e.getMessage());}
        }  
      }      
    } else if (actionCommand.equalsIgnoreCase(ExportPanel.EXPORT_EXCEL)||
        actionCommand.equalsIgnoreCase(ExportPanel.EXPORT_TEXT)){
      ExportOptionsVO expOptions = new ExportOptionsVO(ExportOptionsVO.EXPORT_NO_DEVIATION,null,true,false,false,6);
      if (exportSettings_ != null)
        expOptions = exportSettings_.getSettings();
      if (actionCommand.equalsIgnoreCase(ExportPanel.EXPORT_EXCEL)||actionCommand.equalsIgnoreCase(ExportPanel.EXPORT_TEXT)){
        Hashtable<String,String> expIdToString = new Hashtable<String,String>();
        for (String expId : experimentNames_)
          expIdToString.put(expId, heatMapListener_.getDisplayName(expId));
        if (actionCommand.equalsIgnoreCase(ExportPanel.EXPORT_EXCEL)){
          exportFileChooser_.setFileFilter(new FileNameExtensionFilter("Microsoft Office Excel Woorkbook (*.xlsx)","xlsx"));
          int returnVal = exportFileChooser_.showSaveDialog(new JFrame());
          if (returnVal == JFileChooser.APPROVE_OPTION) {
            File fileToStore = exportFileChooser_.getSelectedFile();
            @SuppressWarnings("rawtypes")
            Vector results = this.checkFileStorage(fileToStore,"xlsx");
            fileToStore = (File)results.get(0);
            if ((Boolean)results.get(1)){
              try {
                int maxIsotope = Integer.parseInt((String)maxIsotopes_.getSelectedItem());
                String preferredUnit = HeatMapDrawing.extractPreferredUnitForExp(preferredUnit_);
                Hashtable<String,Hashtable<String,Vector<Double>>> resultValues = HeatMapDrawing.extractValuesOfInterest(resultsOfOneGroup_, maxIsotope, settingsVO_, preferredUnit, expOptions,modifications_);
//              if (preferredUnit!=null&&preferredUnit.length()>0)
                preferredUnit = StaticUtils.getCorrespondingUnit(settingsVO_,preferredUnit,true);
                BarChartPainter.exportToFile(molGroupName_,new BufferedOutputStream(new FileOutputStream(fileToStore)),true, maxIsotope, getSelectedMoleculeNames(), 
                    experimentNames_,expIdToString, resultValues, preferredUnit, StaticUtils.getAreaTypeString(settingsVO_),expOptions,modifications_);
              }
              catch (NumberFormatException e) {new WarningMessage(new JFrame(), "Error", e.getMessage());}
              catch (FileNotFoundException e) {new WarningMessage(new JFrame(), "Error", e.getMessage());}
              catch (IOException e) {new WarningMessage(new JFrame(), "Error", e.getMessage());}
              catch (CalculationNotPossibleException e) {new WarningMessage(new JFrame(), "Error", e.getMessage());}
            }  
          }
        } else if (actionCommand.equalsIgnoreCase(ExportPanel.EXPORT_TEXT)){
          exportFileChooser_.setFileFilter(new FileNameExtensionFilter("Text (*.txt)","txt"));
          int returnVal = exportFileChooser_.showSaveDialog(new JFrame());
          if (returnVal == JFileChooser.APPROVE_OPTION) {
            File fileToStore = exportFileChooser_.getSelectedFile();
            @SuppressWarnings("rawtypes")
            Vector results = this.checkFileStorage(fileToStore,"txt");
            fileToStore = (File)results.get(0);
            if ((Boolean)results.get(1)){
              try {
                int maxIsotope = Integer.parseInt((String)maxIsotopes_.getSelectedItem());
                String preferredUnit = HeatMapDrawing.extractPreferredUnitForExp(preferredUnit_);
                Hashtable<String,Hashtable<String,Vector<Double>>> resultValues = HeatMapDrawing.extractValuesOfInterest(resultsOfOneGroup_, maxIsotope, settingsVO_, preferredUnit, expOptions,modifications_);
                //if (preferredUnit!=null&&preferredUnit.length()>0)
                preferredUnit = StaticUtils.getCorrespondingUnit(settingsVO_,preferredUnit,true);
                BarChartPainter.exportToFile(molGroupName_,new BufferedOutputStream(new FileOutputStream(fileToStore)),false, maxIsotope, getSelectedMoleculeNames(), 
                    experimentNames_,expIdToString, resultValues, preferredUnit, StaticUtils.getAreaTypeString(settingsVO_),expOptions,modifications_);
              }
              catch (NumberFormatException e) {new WarningMessage(new JFrame(), "Error", e.getMessage());}
              catch (FileNotFoundException e) {new WarningMessage(new JFrame(), "Error", e.getMessage());}
              catch (IOException e) {new WarningMessage(new JFrame(), "Error", e.getMessage());}
              catch (CalculationNotPossibleException e) {new WarningMessage(new JFrame(), "Error", e.getMessage());}
            }  
          }      
        }
      }
    } else if (actionCommand.equalsIgnoreCase(ExportPanel.EXPORT_MZTAB) || actionCommand.equalsIgnoreCase(ExportPanel.EXPORT_RDB)) {
    	// get the directory where to store the mzTab files (1 per experiment)
      exportFileChooser_.setFileSelectionMode(JFileChooser.FILES_ONLY);
      int equalPartsAtBeginning = StaticUtils.detectNameUnequalitiesBeforeAndAfter(experimentNames_)[0];
      String fileName = "";
      String extension = "";
      String confirmDialogTitle = "";
      FileNameExtensionFilter filter = null;
      if (actionCommand.equalsIgnoreCase(ExportPanel.EXPORT_MZTAB)){
        fileName = "export-mztab.txt";
        extension = "mztab.txt";
        filter = new FileNameExtensionFilter("Text (*.txt)","txt");
        confirmDialogTitle = "mzTab export";
      }else if (actionCommand.equalsIgnoreCase(ExportPanel.EXPORT_RDB)){
        fileName = "unified.tab";
        extension = "unified.tab";
        filter = new FileNameExtensionFilter("Tab (*.tab)","tab");
        confirmDialogTitle = "RDB export";
      }
      if (equalPartsAtBeginning>0){
        fileName = experimentNames_.get(0).substring(0,equalPartsAtBeginning);
        if (!fileName.endsWith("-")&&!fileName.endsWith("_")) fileName += "-";
        fileName += extension;
      }
      exportFileChooser_.setSelectedFile(new File(fileName));
      exportFileChooser_.setFileFilter(filter);
      if (JOptionPane.showConfirmDialog(this, "Only selected analytes with the current isotopes selected will be exported! Continue?",confirmDialogTitle,JOptionPane.YES_NO_OPTION) == JOptionPane.YES_OPTION){
        int returnVal = exportFileChooser_.showSaveDialog(new JFrame());	
    	    if (returnVal != JFileChooser.APPROVE_OPTION)
    	      return;
    	
    	    File exportFile = exportFileChooser_.getSelectedFile();
    	    if (actionCommand.equalsIgnoreCase(ExportPanel.EXPORT_MZTAB)) heatMapListener_.exportMzTab(exportFile);
    	    else if (actionCommand.equalsIgnoreCase(ExportPanel.EXPORT_RDB)) heatMapListener_.exportRdb(exportFile);
      }
      exportFileChooser_.setSelectedFile(new File(""));
    	
    } else if (actionCommand.equalsIgnoreCase(ExportPanel.EXPORT_CHROMS)){
      chromExport_.setVisible(true);
    } else if (actionCommand.equalsIgnoreCase(ExportPanel.EXPORT_MAF)) {
    // get the directory where to store the maf.tsv file (1 per experiment)
      exportFileChooser_.setFileSelectionMode(JFileChooser.FILES_ONLY);
      int equalPartsAtBeginning = StaticUtils.detectNameUnequalitiesBeforeAndAfter(experimentNames_)[0];
      String fileName = "export-maf.tsv";
      if (equalPartsAtBeginning>0){
        fileName = experimentNames_.get(0).substring(0,equalPartsAtBeginning);
        if (!fileName.endsWith("-")&&!fileName.endsWith("_")) fileName += "-";
        fileName += "maf.tsv";
      }
      exportFileChooser_.setSelectedFile(new File(fileName));
      exportFileChooser_.setFileFilter(new FileNameExtensionFilter("Tab-separated (*.tsv)","tsv"));
      if (JOptionPane.showConfirmDialog(this, "Only selected analytes with the current isotopes selected will be exported! Continue?","mzTab export",JOptionPane.YES_NO_OPTION) == JOptionPane.YES_OPTION){
        int returnVal = exportFileChooser_.showSaveDialog(new JFrame());
        if (returnVal != JFileChooser.APPROVE_OPTION)
          return;    
        File exportFile = exportFileChooser_.getSelectedFile();
        heatMapListener_.exportMaf(exportFile);
      }
      exportFileChooser_.setSelectedFile(new File(""));
    

    } else if (actionCommand.equalsIgnoreCase("AcceptChromExport")){
      Vector<String> expsToExport = new Vector<String>();
      Vector<String> analsToExport = new Vector<String>();
      Vector<String> resultFiles = new Vector<String>();
      Vector<String> chromsToUse = new Vector<String>();
      Vector<String> modsToUse = new Vector<String>();
      for (int i=0;i!=this.experimentNames_.size();i++ ){
        String name = this.experimentNames_.get(i);
        if (chromExport_.isExperimentSelected(heatMapListener_.getDisplayName(name))){
          if (compVOs_[i]!=null && compVOs_[i].length>0){
            String resultFile = compVOs_[i][0].getAbsoluteFilePath();
            resultFiles.add(resultFile);
            String chromFileBase = StaticUtils.extractChromBaseName(resultFile,name);
            boolean chromFileExists = false;
            if (chromFileBase!=null && chromFileBase.length()>0){
              chromFileExists = true;        
            }
            if (chromFileExists && !StaticUtils.existChromFiles(chromFileBase)) chromFileExists = false;           
            if (chromFileExists){
              String chromFilePath = chromFileBase+".chrom";
              chromsToUse.add(chromFilePath);
              expsToExport.add(heatMapListener_.getDisplayName(name));
            }
          }  
        }
      }
      for (String mod : this.modifications_){
        if (chromExport_.isModificationSelected(mod))modsToUse.add(mod);
      }
      for (String name : this.moleculeNames_){
        if (chromExport_.isMoleculeSelected(name)){
          if (modifications_.size()==1)
            analsToExport.add(name);
          else if (modifications_.size()>1){
            for (String mod : modsToUse) analsToExport.add(name+"_"+mod);
          }
        }
        
      }
      if (expsToExport.size()>0 && analsToExport.size()>0 && chromsToUse.size()>0 && modsToUse.size()>0){
        exportFileChooser_.setSelectedFile(null);
        exportFileChooser_.updateUI();
        if (chromExport_.getPictureType().equalsIgnoreCase(ExportPanel.EXPORT_PNG))
          exportFileChooser_.setFileFilter(new FileNameExtensionFilter("PNG (*.png)","png"));
        else if (chromExport_.getPictureType().equalsIgnoreCase(ExportPanel.EXPORT_SVG))
          exportFileChooser_.setFileFilter(new FileNameExtensionFilter("SVG (*.svg)","svg"));
        int returnVal = exportFileChooser_.showSaveDialog(new JFrame());
        if (returnVal == JFileChooser.APPROVE_OPTION) {
          File fileToStore = exportFileChooser_.getSelectedFile();
          @SuppressWarnings("rawtypes")
          Vector results = null;
          if (chromExport_.getPictureType().equalsIgnoreCase(ExportPanel.EXPORT_PNG))
            results = this.checkFileStorage(fileToStore,"png");
          if (chromExport_.getPictureType().equalsIgnoreCase(ExportPanel.EXPORT_SVG))
            results = this.checkFileStorage(fileToStore,"svg");
          fileToStore = (File)results.get(0);
          if ((Boolean)results.get(1) && chromExportThread_==null){
            chromExportThread_ = new ChromExportThread(groupName_,fileToStore,chromExport_.getPictureDimension(),
                expsToExport,resultFiles,chromsToUse,analsToExport,chromExport_.getPictureType(),exportSettings_.getSettings().isAnalyteInColumn(),
                rtTolerance_,isLookup_,esLookup_);
            chromExportThread_.start();
            spinnerLabel_.setVisible(true);
            exportProgressPanel_.setVisible(true);
          } else if (chromExportThread_!=null)
            new WarningMessage(new JFrame(), "Error", "There is already a chromatogram export running! Wait until it is finished!");
        }
      }
    } else if (actionCommand.equalsIgnoreCase("Choose just one peak for doubles") || actionCommand.equalsIgnoreCase("Quant. anal. at not found")||
        actionCommand.equalsIgnoreCase("Take exact peak for others")){
      String expName = experimentNames_.get(lastClickedCellPos_[0]);
      String molName = moleculeNames_.get(lastClickedCellPos_[1]);
      ResultCompVO vo = compVOs_[lastClickedCellPos_[0]][lastClickedCellPos_[1]];
      int maxIsotopes = Integer.parseInt((String)maxIsotopes_.getSelectedItem());
      Vector<String> foundUpdateables = new Vector<String>();
//      Hashtable<String,String> updateableAndAnalyteBefore = new Hashtable<String,String>();
      Vector<AutoAnalyteAddVO> updateableAndAnalyteBefore = new Vector<AutoAnalyteAddVO>();
      Hashtable<String,String> modHash = new Hashtable<String,String>();
      int maxIsotope = 0;
      for (int i=0;i!=experimentNames_.size();i++){
        ResultCompVO otherVO = compVOs_[i][lastClickedCellPos_[1]];
        if (actionCommand.equalsIgnoreCase("Choose just one peak for doubles")){
          if (i!=lastClickedCellPos_[0] && otherVO.getMoreThanOnePeak(otherVO.getAvailableIsotopeNr(maxIsotopes))){
            foundUpdateables.add(otherVO.getAbsoluteFilePath());
            for (String modName : modifications_){
              if (otherVO.getMoreThanOnePeak(otherVO.getAvailableIsotopeNr(maxIsotopes), modName)) modHash.put(modName, modName);
            }
          }
        } else if (actionCommand.equalsIgnoreCase("Quant. anal. at not found")||actionCommand.equalsIgnoreCase("Take exact peak for others")){
          if (i!=lastClickedCellPos_[0] && !otherVO.existsInFile()){
            String previousElement = "";
            for (int j=lastClickedCellPos_[1]-1;j>-1;j--){
              if (compVOs_[i][j].existsInFile()){
                previousElement = moleculeNames_.get(j);
                break;
              }
            }
            updateableAndAnalyteBefore.add(new AutoAnalyteAddVO(experimentNames_.get(i),otherVO.getAbsoluteFilePath(), previousElement));
//            updateableAndAnalyteBefore.put(otherVO.getAbsoluteFilePath(), previousElement);
            maxIsotope = vo.getUsedIsotpes();
            for (String modName : modifications_){
              if (vo.containsMod(modName))modHash.put(modName, modName);
            }
          }
        }
      }
      if (foundUpdateables.size()>0 || updateableAndAnalyteBefore.size()>0){
        String confirmationMessage = "Are you sure you want to eliminate double peaks for "+molName+"? The peak nearest to the RT of "+expName+" is taken!";
        if (actionCommand.equalsIgnoreCase("Quant. anal. at not found")){
          confirmationMessage = "Are you sure you want to automatically add peaks for "+molName+"? The peak nearest to the RT of "+expName+" is taken!";
        }else if (actionCommand.equalsIgnoreCase("Take exact peak for others")){
          confirmationMessage = "Are you sure you want to automatically add peaks for "+molName+"? The exact peak borders of "+expName+" are taken!";
        }
        Vector<String> selectedMods = CheckBoxOptionPane.showConfirmDialog(new JFrame(), "Confirmation", confirmationMessage, new Vector<String>(modHash.values()),true);
        if (selectedMods.size()>0){
          boolean exactProbePosition = false;
          if (actionCommand.equalsIgnoreCase("Take exact peak for others")) exactProbePosition = true;
          if (actionCommand.equalsIgnoreCase("Choose just one peak for doubles"))
            heatMapListener_.eliminateDoublePeaks(groupName_,molName,vo.getAbsoluteFilePath(),selectedMods,foundUpdateables);
          else if (actionCommand.equalsIgnoreCase("Quant. anal. at not found")||actionCommand.equalsIgnoreCase("Take exact peak for others"))
            heatMapListener_.addAnalyteEverywhereAtPosition(groupName_,molName,vo.getAbsoluteFilePath(),selectedMods,updateableAndAnalyteBefore,maxIsotope,exactProbePosition);          
        }
      } else {
        if (actionCommand.equalsIgnoreCase("Choose just one peak for doubles"))
          new WarningMessage(new JFrame(), "Error", "There are no double peaks to eliminate for "+molName+"!");
        else if (actionCommand.equalsIgnoreCase("Quant. anal. at not found"))
          new WarningMessage(new JFrame(), "Error", "There are no  peaks to add for "+molName+"!");
      }
    } else if (actionCommand.equalsIgnoreCase("Remove analyte in all probes")){
      Hashtable<String,String> selectedAnalytes = new Hashtable<String,String>(selectedMolecules_);
      selectedAnalytes.put(lastClickedAnalyte_, lastClickedAnalyte_);
      String analyteList = "";
      for (String analName : selectedAnalytes.keySet()) analyteList += analName+", ";
      Vector<String> selectedMods = CheckBoxOptionPane.showConfirmDialog(new JFrame(), "Confirmation", "Do you really want to delete "+groupName_+" "+analyteList+" in all of the probes?", modifications_,false);
      if (selectedMods.size()>0){
        Hashtable<String,String> foundUpdateables = new Hashtable<String,String>();
        Hashtable<Integer,Integer> rowsWhereFound = new Hashtable<Integer,Integer>();
        for (int i=0;i!=this.moleculeNames_.size();i++){
          if (selectedAnalytes.containsKey(moleculeNames_.get(i)))
            rowsWhereFound.put(i, i);
        }
        
        for (Integer rowWhereFound: rowsWhereFound.keySet()){
          for (int i=0;i!=experimentNames_.size();i++){
            ResultCompVO otherVO = compVOs_[i][rowWhereFound];
            if (otherVO.getAbsoluteFilePath()!=null && otherVO.getAbsoluteFilePath().length()>0 && otherVO.existsInFile()){
              boolean foundMod = false;
              for (String modName : selectedMods){
                if (otherVO.containsMod(modName))foundMod = true; 
              }
              if (foundMod)
                foundUpdateables.put(otherVO.getAbsoluteFilePath(),otherVO.getAbsoluteFilePath());
            }  
          }
        }
        heatMapListener_.eliminateAnalyteEverywhere(groupName_, selectedAnalytes, selectedMods, new Vector<String>(foundUpdateables.values()));
      }
    }else if (actionCommand.equalsIgnoreCase("Select analyte") || actionCommand.equalsIgnoreCase("Deselect analyte")){  
      Graphics2D g2 = (Graphics2D)renderedImage_.getGraphics();
      this.selectAnalyte(lastClickedAnalyte_);
      if (actionCommand.equalsIgnoreCase("Select analyte")){
        this.selectAnalyte(lastClickedAnalyte_);
        g2.setColor(Color.BLUE);
      }else if (actionCommand.equalsIgnoreCase("Deselect analyte")){
        g2.setColor(Color.BLACK);
        this.deselectAnalyte(lastClickedAnalyte_);
      }  
      Rectangle rectForName = heatmap_.getRectangleForRowName(lastClickedAnalyte_);
      rectForName.setLocation(new Point(rectForName.getLocation().x+imagePositionX_,rectForName.getLocation().y+imagePositionY_));
      g2.fillRect(rectForName.x+1,rectForName.y+1, rectForName.width-2, rectForName.height-2);
      Font descriptionFont = new Font("Dialog",Font.PLAIN, 9);
      FontMetrics descriptionFontMetrics = g2.getFontMetrics();
      int textHeight = descriptionFontMetrics.getHeight();
      g2.setColor(Color.WHITE);
      g2.setFont(descriptionFont);
      g2.drawString(lastClickedAnalyte_,rectForName.x+1,rectForName.y + textHeight/2);
      this.invalidate();
      this.updateUI();
    } else if (actionCommand.equalsIgnoreCase("stopChromExport")){
      chromExportThread_.stopThread();
    }
  }
  
  private static Hashtable<String,Hashtable<String,Vector<Double>>> extractValuesOfInterest(Hashtable<String,Hashtable<String,ResultCompVO>> vos, int maxIsotope, ResultDisplaySettingsVO settingVO, String preferredUnit, ExportOptionsVO expOptions, Vector<String> modifications) throws CalculationNotPossibleException{
    Hashtable<String,Hashtable<String,Vector<Double>>> results = new Hashtable<String,Hashtable<String,Vector<Double>>>();
    for (String molKey : vos.keySet()){
      Hashtable<String,Vector<Double>> expValues = new Hashtable<String,Vector<Double>>();
      for (String expKey : vos.get(molKey).keySet()){
        ResultCompVO compVO = vos.get(molKey).get(expKey);
        int isoNr = maxIsotope;
        if (!settingVO.getType().equalsIgnoreCase("relative to measured class amount") && 
            !settingVO.getType().equalsIgnoreCase("relative to total amount"))
          isoNr = compVO.getAvailableIsotopeNr(maxIsotope);
        double myArea = compVO.getArea(isoNr, settingVO);
        if (settingVO.getType().equalsIgnoreCase("relative to measured class amount") || 
            settingVO.getType().equalsIgnoreCase("relative to total amount"))
          compVO.getArea(maxIsotope, settingVO);
        myArea = StaticUtils.getAreaInCorrespondingUnit(myArea, preferredUnit);
        Vector<Double> areaPlusDev = new Vector<Double>();
        areaPlusDev.add(myArea);
        if (expOptions!=null&&expOptions.getExportType()!=ExportOptionsVO.EXPORT_NO_DEVIATION){
          ResultCompGroupVO groupVO = (ResultCompGroupVO)compVO;
          double sdValue = -1d;
          if (expOptions.getExportType() == ExportOptionsVO.EXPORT_SD_DEVIATION || expOptions.getExportType() == ExportOptionsVO.EXPORT_SD_DEV_AND_ERROR){
            sdValue = groupVO.getAreaSD(isoNr, settingVO)*Double.parseDouble(expOptions.getSdValue());
            areaPlusDev.add(getExportableSdValue(sdValue,preferredUnit));
          }
          if (expOptions.getExportType() == ExportOptionsVO.EXPORT_SD_ERROR || expOptions.getExportType() == ExportOptionsVO.EXPORT_SD_DEV_AND_ERROR){
            sdValue = groupVO.getAreaSE(isoNr, settingVO);
            areaPlusDev.add(getExportableSdValue(sdValue,preferredUnit));
          }
        }
        if (expOptions!=null&&expOptions.isExportRT()){
          for (String modName : modifications){
            areaPlusDev.add(compVO.getRetentionTime(modName));
            if (expOptions!=null&&expOptions.isExportRTDev()){
              ResultCompGroupVO groupVO = (ResultCompGroupVO)compVO;
              areaPlusDev.add(groupVO.getRetentionTimeSD(modName));
            }
          }
        }
        expValues.put(expKey, areaPlusDev);
      }   
      results.put(molKey, expValues);
    }
    return results;
  }
  
  private static double getExportableSdValue(double sdValue, String preferredUnit){
    double sd = -1d;
    if (Double.isInfinite(sdValue) || Double.isNaN(sdValue))
      sd = -1d;
    else
      sd = StaticUtils.getAreaInCorrespondingUnit(sdValue, preferredUnit);
    return sd;
  }
  
  @SuppressWarnings({ "unchecked", "rawtypes" })
  private Vector checkFileStorage (File file, String suffix){
    Vector results = new Vector();
    boolean store = true;
    File fileToStore = new File(file.getAbsolutePath());
    if (fileToStore.exists()){
      if (JOptionPane.showConfirmDialog(this, "The file "+fileToStore.getName()+" exists! Replace existing file?") != JOptionPane.YES_OPTION)
        store = false;
    }else{
      if (fileToStore.getName().indexOf(".")==-1){
        fileToStore = new File(fileToStore.getAbsoluteFile()+"."+suffix);
        if (fileToStore.exists()){
          if (JOptionPane.showConfirmDialog(this, "The file "+fileToStore.getName()+" exists! Replace existing file?") != JOptionPane.YES_OPTION)
            store = false;
        }  
      }  
    }
    results.add(fileToStore);
    results.add(store);
    return results;
  }
  
  private class MyListener extends MouseInputAdapter {
    public void mouseMoved(MouseEvent e) {
      int x = e.getX();
      int y = e.getY();
      if (heatmap_!=null && renderedImage_!=null) {
        Graphics2D g2 = (Graphics2D)e.getComponent().getGraphics();
        if (rectToDraw_!=null){
          g2.setColor(Color.BLACK);
          g2.drawRect(rectToDraw_.x,rectToDraw_.y, rectToDraw_.width, rectToDraw_.height);
          if (statusText_!=null){
            statusText_.setText("");
          }  
        }
        if ((y>=imagePositionY_) && (y<imagePositionY_+renderedImage_.getHeight()) && (x>=imagePositionX_) && (x<(imagePositionX_+renderedImage_.getWidth()))) {
          int xInImage = x-imagePositionX_;
          int yInImage = y-imagePositionY_;
          if (heatmap_.getExpressionImageXStart()<=xInImage&&xInImage<heatmap_.getExpressionImageXEnd()&&
              heatmap_.getExpressionImageYStart()<=yInImage&&yInImage<heatmap_.getExpressionImageYEnd()){
            Vector<String> cellName = heatmap_.getCellName(xInImage, yInImage);
            if (cellName!=null&&cellName.size()==2&&cellName.get(0)!=null&&cellName.get(1)!=null){
//              System.out.println(cellName.get(0)+" ; "+cellName.get(1));
              int[] cellPos = heatmap_.getCellPosition(xInImage, yInImage);
              if (statusText_!=null){
                ResultCompVO compVO = compVOs_[cellPos[0]][cellPos[1]];
                int maxIsotope = Integer.parseInt((String)maxIsotopes_.getSelectedItem());
                String statusText = "";
                try {
                  double relativeValue = compVO.getRelativeValue(compVO.getAvailableIsotopeNr(maxIsotope),settingsVO_);
                  String sampleTypeText = "Sample";
                  if (isGrouped_) sampleTypeText = "Group";
                  statusText = "Lipid: "+cellName.get(1)+"; "+sampleTypeText+": "+cellName.get(0)+"; Relative: "+StaticUtils.extractDisplayValue(relativeValue);
                  double value = compVO.getArea(compVO.getAvailableIsotopeNr(maxIsotope), settingsVO_);
                  if (settingsVO_.getType().equalsIgnoreCase("relative to measured class amount") || 
                      settingsVO_.getType().equalsIgnoreCase("relative to total amount"))
                    value = compVO.getArea(maxIsotope, settingsVO_);
                  String selectedAreaString = StaticUtils.getAreaString(value, settingsVO_, preferredUnit_.get(cellName.get(1)));
                  if (selectedAreaString!=null&&selectedAreaString.length()>0)
                    statusText += "; "+selectedAreaString;
                  double standValue = compVO.getStandardizedArea(compVO.getAvailableIsotopeNr(maxIsotope), settingsVO_.getISStandMethod(),  settingsVO_.getESStandMethod(), settingsVO_.considerDilution());
                  statusText += "; Stand-Area: "+standValue;
                  statusText += "; Quant-Area: "+compVO.getOriginalArea(compVO.getAvailableIsotopeNr(maxIsotope))+"; Isos: "+(compVO.getAvailableIsotopeNr(maxIsotope)+1);
                }
                catch (CalculationNotPossibleException e1) {
                  // TODO Auto-generated catch block
                  e1.printStackTrace();
                }
                
                statusText_.setText(statusText);
              }  
              rectToDraw_ = heatmap_.getRectangleForCell(xInImage, yInImage);
              this.drawARectangle(g2);
            }
          } else if (heatmap_.getRowNameStart()<xInImage&& xInImage<heatmap_.getRowNameEnd() &&
              heatmap_.getExpressionImageYStart()<=yInImage&&yInImage<heatmap_.getExpressionImageYEnd()){
            String analyteName = heatmap_.getRowName(x,y);
//            if (analyteName!=null&&analyteName.length()>0&&!isAnalyteSelected(analyteName)){
              rectToDraw_ = heatmap_.getRectangleForRowName(analyteName);
              this.drawARectangle(g2);
//            }else
//              rectToDraw_ = null;
          } else if (heatmap_.getExpressionImageXStart()<=xInImage&&xInImage<heatmap_.getExpressionImageXEnd()&&
              heatmap_.getColumnNameStart()<=yInImage&&yInImage<heatmap_.getColumnNameEnd()){
            String expName = heatmap_.getColumnName(x,y);
            if (expName!=null&&expName.length()>0){
              rectToDraw_ = heatmap_.getRectangleForColumnName(expName);
              this.drawARectangle(g2);              
            }
          }
        }
      }  
    }
    
    private void drawARectangle(Graphics2D g2){
      rectToDraw_.setLocation(new Point(rectToDraw_.getLocation().x+imagePositionX_,rectToDraw_.getLocation().y+imagePositionY_));
      g2.setColor(Color.WHITE);
      g2.drawRect(rectToDraw_.x,rectToDraw_.y, rectToDraw_.width, rectToDraw_.height);      
    }
    
    public void mouseClicked(MouseEvent e) {
      int x = e.getX();
      int y = e.getY();
      if (heatmap_!=null && renderedImage_!=null) {
        Graphics2D g2 = (Graphics2D)e.getComponent().getGraphics();
        if (rectToDraw_!=null){
          g2.setColor(Color.BLACK);
          g2.drawRect(rectToDraw_.x,rectToDraw_.y, rectToDraw_.width, rectToDraw_.height);
          if (statusText_!=null){
            statusText_.setText("");
          }  
        }
        int maxIsotopes = Integer.parseInt((String)maxIsotopes_.getSelectedItem());
        if ((y>=imagePositionY_) && (y<imagePositionY_+renderedImage_.getHeight()) && (x>=imagePositionX_) && (x<(imagePositionX_+renderedImage_.getWidth()))) {
          int xInImage = x-imagePositionX_;
          int yInImage = y-imagePositionY_;
          if (heatmap_.getExpressionImageXStart()<=xInImage&&xInImage<heatmap_.getExpressionImageXEnd()&&
              heatmap_.getExpressionImageYStart()<=yInImage&&yInImage<heatmap_.getExpressionImageYEnd()){
            Vector<String> cellName = heatmap_.getCellName(xInImage, yInImage);
            if (cellName!=null&&cellName.size()==2&&cellName.get(0)!=null&&cellName.get(1)!=null){
              int[] cellPos = heatmap_.getCellPosition(xInImage, yInImage);
              ResultCompVO compVO = compVOs_[cellPos[0]][cellPos[1]];
              if (compVO.getOriginalArea(0)>0 && heatMapListener_!=null){
                if (!isGrouped_){
                  if (e.getButton()==MouseEvent.BUTTON1){
                    if (!heatMapListener_.heatMapClicked(cellName.get(0),compVO.getAbsoluteFilePath(),cellName.get(1)))
                      this.mouseMoved(e);
                  }else if (e.getButton()==MouseEvent.BUTTON3 || e.isPopupTrigger()){
                    if (attentionProbes_.containsKey(cellName.get(0))&&attentionProbes_.get(cellName.get(0)).containsKey(cellName.get(1)))
                      this.mouseMoved(e);
                    else{
                      applySettingsPopup_.show(e.getComponent(), e.getX(), e.getY());
                      lastClickedCellPos_ = cellPos;
//                      lastClickedExp_ = cellName.get(0);
//                      lastClickedMol_ = cellName.get(1);
//                      lastClickedResultVO_ = compVO;
                    }  
                  }
                }else
                  this.mouseMoved(e);
              }else
                this.mouseMoved(e);
            }
          }else if (heatmap_.getRowNameStart()<xInImage&& xInImage<heatmap_.getRowNameEnd() &&
              heatmap_.getExpressionImageYStart()<=yInImage&&yInImage<heatmap_.getExpressionImageYEnd()){
            
            String analyteName = heatmap_.getRowName(x,y);
            boolean returnValue = false;
            if (isGrouped_)
              returnValue = heatMapListener_.analyteGroupClicked(analyteName,groupName_,maxIsotopes,settingsVO_,preferredUnit_.get(analyteName),StaticUtils.getCorrespondingUnit(settingsVO_, preferredUnit_.get(analyteName),true));
            else{
              if (SwingUtilities.isLeftMouseButton(e)){
                returnValue = heatMapListener_.analyteClicked(analyteName,groupName_,maxIsotopes,settingsVO_,preferredUnit_.get(analyteName),StaticUtils.getCorrespondingUnit(settingsVO_, preferredUnit_.get(analyteName),true));
              } else if (SwingUtilities.isRightMouseButton(e)){
                if (isAnalyteSelected(analyteName)){
                  selectItem_.setEnabled(false);
                  deselectItem_.setEnabled(true);
                }else{
                  selectItem_.setEnabled(true);
                  deselectItem_.setEnabled(false);                  
                }                  
                removeAnalytePopup_.show(e.getComponent(), e.getX(), e.getY());
                lastClickedAnalyte_ = analyteName;
              }
            } 
            if (!returnValue)
              this.mouseMoved(e);
          } else if (heatmap_.getExpressionImageXStart()<=xInImage&&xInImage<heatmap_.getExpressionImageXEnd()&&
              heatmap_.getColumnNameStart()<=yInImage&&yInImage<heatmap_.getColumnNameEnd()){
            String expName = heatmap_.getColumnName(x,y);
            if (e.getButton()==MouseEvent.BUTTON1){
              boolean returnValue = false;
              String preferredUnit = extractPreferredUnitForExp(preferredUnit_);
              if (isGrouped_)
                returnValue = heatMapListener_.experimentGroupClicked(expName,groupName_,maxIsotopes,settingsVO_,preferredUnit,StaticUtils.getCorrespondingUnit(settingsVO_, preferredUnit,true));
              else
                returnValue = heatMapListener_.experimentClicked(expName,groupName_,maxIsotopes,settingsVO_,preferredUnit,StaticUtils.getCorrespondingUnit(settingsVO_, preferredUnit,true));
              if (!returnValue)
                this.mouseMoved(e);
            }else if (!isGrouped_ && e.getButton()==MouseEvent.BUTTON3){
              InputDialog dlg = new InputDialog(new JFrame(), "", "Sample name",heatMapListener_.getDisplayName(expName));
              heatMapListener_.setDisplayName(expName, dlg.getEnteredText());
            }  
          }
        }
      }        
    }
  }
  
  private Vector<String> getUnselectedMoleculeNames(){
    Vector<String> unselectedNames = new Vector<String>();
    for (String name : moleculeNames_){
      if (!this.selectionSettings_.isSelected(name)) unselectedNames.add(name);
    }
    return unselectedNames;    
  }
  
  public Vector<String> getSelectedMoleculeNames(){
    
    Vector<String> selectedNames = new Vector<String>();
    for (String name : moleculeNames_){
      boolean accept = true;
      if (this.selectionSettings_.isSelected(name)){
        if (!this.showInternalStandards_.isSelected() && isLookup_.containsKey(name))
          accept = false;
        else if (!this.showExternalStandards_.isSelected() && esLookup_.containsKey(name))
          accept = false;
      }else
        accept = false;
      if (accept)
        selectedNames.add(name);
    }
    return selectedNames;
  }

//  public boolean hasNoInternalStandards(){
//    return !(this.showInternalStandards_.isSelected());
//  }
  
  private class SelectionItemListener implements java.awt.event.ItemListener
  {
    private String m_ctrl;
    
    public SelectionItemListener(String ctrl){
      m_ctrl = ctrl;
    }

    public void itemStateChanged(ItemEvent e)  {
      if (m_ctrl.equalsIgnoreCase("ChangeIsotope")){
        if (e.getStateChange()==ItemEvent.SELECTED){
          if (parentAction_);
            heatMapListener_.changeIsotopesUsed(groupName_, isGrouped_, Integer.parseInt((String)maxIsotopes_.getSelectedItem()));
          generateHeatMap();
          invalidate();
          updateUI();
        }
      }
    }  
  }
  

  
  

  

  
  private static String extractPreferredUnitForExp(Hashtable<String,String> preferredUnits){
    String preferredUnit = "";
    List<Integer> unitValues = new ArrayList<Integer>();
    for (String unit : preferredUnits.values()){
      if (unit.equalsIgnoreCase(""))
        unitValues.add(0);
      else if (unit.equalsIgnoreCase("%"))
        unitValues.add(1);
      else if (unit.equalsIgnoreCase("\u2030"))
        unitValues.add(2);      
      else if (unit.equalsIgnoreCase("m"))
        unitValues.add(3);
      else if (unit.equalsIgnoreCase("\u03BC"))
        unitValues.add(4);
      else if (unit.equalsIgnoreCase("n"))
        unitValues.add(5);
      else if (unit.equalsIgnoreCase("p"))
        unitValues.add(6);
      else if (unit.equalsIgnoreCase("f"))
        unitValues.add(7);
      else if (unit.equalsIgnoreCase("a"))
        unitValues.add(8);
    }
    Collections.sort(unitValues);
    int unitPlaceHolder = 0;
    if (unitValues.size()%2 == 0){
      unitPlaceHolder = unitValues.get(((unitValues.size()) / 2));
    }else{
      unitPlaceHolder = unitValues.get(((unitValues.size() + 1) / 2) - 1);
    }
    if (unitPlaceHolder == 0)
      preferredUnit = "";
    if (unitPlaceHolder == 1)
      preferredUnit = "%";
    if (unitPlaceHolder == 2)
      preferredUnit = "\u2030";
    if (unitPlaceHolder == 3)
      preferredUnit = "m";
    if (unitPlaceHolder == 4)
      preferredUnit = "\u03BC";
    if (unitPlaceHolder == 5)
        preferredUnit = "n";
    if (unitPlaceHolder == 6)
      preferredUnit = "p";
    if (unitPlaceHolder == 7)
      preferredUnit = "f";
    if (unitPlaceHolder == 8)
      preferredUnit = "a";
    return preferredUnit;
  }
  
  public void setISSelected(boolean selected){
    this.showInternalStandards_.setSelected(selected);
    parentAction_ = false;
    this.actionPerformed(CHANGE_IS_STATUS);
    parentAction_ = true;
  }
  
  public void setESSelected(boolean selected){
    this.showExternalStandards_.setSelected(selected);
    parentAction_ = false;
    this.actionPerformed(CHANGE_ES_STATUS);
    parentAction_ = true;
  }
  
  public void setSelectedIsotope(int isotope){
    parentAction_ = false;
    maxIsotopes_.setSelectedItem(String.valueOf(isotope));
    parentAction_ = true;
  }

  public void deselectAnalyte(String analyteName)
  {
    selectedMolecules_.remove(analyteName);
  }

  public boolean isAnalyteSelected(String analyteName)
  {
    return selectedMolecules_.containsKey(analyteName);
  }

  public void selectAnalyte(String analyteName)
  {
    selectedMolecules_.put(analyteName, analyteName);
    
  }
  
  public int getSelectedIsotope(){
    int isotope = 0;
    if (maxIsotopes_.getSelectedItem()!=null)
      isotope = new Integer((String)maxIsotopes_.getSelectedItem()); 
    return isotope;
  }
  
  private Vector<String> getDisplayNames(){
    Vector<String> names = new Vector<String>();
    fromShortToExpName_ = new Hashtable<String,String>();
    for (String name : this.experimentNames_){
      String displayName = heatMapListener_.getDisplayName(name);
      names.add(displayName);
      fromShortToExpName_.put(displayName, name);
    }
    return names;
  }
  
  private Vector<String> getISs(){
    Vector<String> iss = new Vector<String>();
    for (String anal : moleculeNames_){
      if (this.isLookup_.containsKey(anal))iss.add(anal);
    }
    return iss;
  }
  
  private Vector<String> getESs(){
    Vector<String> ess = new Vector<String>();
    for (String anal : moleculeNames_){
      if (this.esLookup_.containsKey(anal))ess.add(anal);
    }
    return ess;
  }
  
  private void initTimer(){
    timer_ = new java.util.Timer();
    timer_.schedule(new ThreadSupervisor(), 10, 1000);
  }
  
  private class ThreadSupervisor extends TimerTask{

    public void run()
    {
      handleTimerEvent();
    }    
  }

  private void handleTimerEvent(){
    if (this.chromExportThread_!=null && !this.chromExportThread_.finished() && 
        (this.chromExportThread_.getErrorString()==null||this.chromExportThread_.getErrorString().length()==0)){
      if (chromExportThread_.getTotalAmountOfLipids()>0&&chromExportThread_.getCurrentLipidCount()>0){
        exportLabel_.setText("Chrom-export: "+chromExportThread_.getCurrentExperiment()+": "+chromExportThread_.getCurrentLipid()+" ("+chromExportThread_.getCurrentLipidCount()+"/"+chromExportThread_.getTotalAmountOfLipids()+")");
        exportProgress_.setValue(((chromExportThread_.getCurrentLipidCount()-1)*100)/chromExportThread_.getTotalAmountOfLipids());
      }
    }
    if (this.chromExportThread_!=null && this.chromExportThread_.finished()){
      if (this.chromExportThread_.getErrorString()!=null&&this.chromExportThread_.getErrorString().length()>0){
        spinnerLabel_.setVisible(false);
        new WarningMessage(new JFrame(), "Error", chromExportThread_.getErrorString());
        exportProgressPanel_.setVisible(false);
      }
      if (this.chromExportThread_.getErrorString()==null||this.chromExportThread_.getErrorString().length()==0){
        exportProgress_.setValue(100);
        exportLabel_.setText("Finished");
        spinnerLabel_.setVisible(false);
        new WarningMessage(new JFrame(), "Information", "The chrom export is finished!");
      }  
      chromExportThread_ = null;
      exportProgressPanel_.setVisible(false);

    }
    
    if (this.chromExportThread_!=null && !this.chromExportThread_.finished() && 
        (this.chromExportThread_.getErrorString()!=null&&this.chromExportThread_.getErrorString().length()>0)){
      spinnerLabel_.setVisible(false);
      new WarningMessage(new JFrame(), "Error", chromExportThread_.getErrorString());
      exportProgressPanel_.setVisible(false);      
    }
  }
  
  public int getPossibleIsotopeNumber(String molecule){
    int maxIsotope = Integer.parseInt((String)maxIsotopes_.getSelectedItem());
    return StaticUtils.getMaxApplicableIsotope(resultsOfOneGroup_.get(molecule),maxIsotope);
  }
  
  public void cleanup(){
    if (heatmap_!=null) heatmap_.cleanup();
    heatmap_ = null;
    if (renderedImage_!=null) renderedImage_.flush();
    renderedImage_ = null;
    rectToDraw_=null;
    if (label_!= null)label_.getGraphics().dispose();
    label_ = null;
    if (heatMapPanel_!=null) heatMapPanel_.getGraphics().dispose();
    heatMapPanel_ = null;
    if (statusText_!=null) statusText_.getGraphics().dispose();
    statusText_ = null;
    if (exportProgressPanel_!=null) exportProgressPanel_.getGraphics().dispose();
    exportProgressPanel_ = null;
    if (compVOs_!=null){
      for (int i=0;i!=compVOs_.length;i++){
        for (int j=0;j!=compVOs_[i].length;j++){
          compVOs_[i][j].cleanup();
          compVOs_[i][j] = null;
        }
      }
    }
    compVOs_ = null;
    heatMapListener_ = null;
    groupName_ = null;
    if (showInternalStandards_!=null) showInternalStandards_.getGraphics().dispose();
    showInternalStandards_ = null;
    if (showInternalStandards_!=null) showInternalStandards_.getGraphics().dispose();    
    showExternalStandards_ = null;
    if (markDoublePeaks_!=null) markDoublePeaks_.getGraphics().dispose();
    markDoublePeaks_ = null;
    if (maxIsotopes_!=null) maxIsotopes_.getGraphics().dispose();
    maxIsotopes_ = null;
    if (resultsOfOneGroup_!=null){
      for (Hashtable<String,ResultCompVO> results : resultsOfOneGroup_.values()){
        for (ResultCompVO vo : results.values()){
          vo.cleanup();
          vo = null;
        }
      }
    }
    resultsOfOneGroup_ = null;
    if (experimentNames_!=null) experimentNames_.clear();
    experimentNames_ = null;
    if (moleculeNames_!=null) moleculeNames_.clear();
    moleculeNames_ = null;
    if (isLookup_!=null) isLookup_.clear();
    isLookup_ = null;
    if (esLookup_!=null) esLookup_.clear();
    esLookup_ = null;
    settingsVO_ = null;
    displaySettings_ = null;
    selectionSettings_ = null;
    exportSettings_ = null;
    combinedChartSettings_ = null;
    if (preferredUnit_!=null) preferredUnit_.clear();
    preferredUnit_ = null;
    if (exportFileChooser_!=null && exportFileChooser_.getGraphics()!=null) exportFileChooser_.getGraphics().dispose();
    exportFileChooser_ = null;
    molGroupName_ = null;
    if (exportLabel_!=null) exportLabel_.getGraphics().dispose();
    exportLabel_ = null;
    if (exportProgress_!=null) exportProgress_.getGraphics().dispose();
    exportProgress_ = null;
    if (spinnerLabel_!=null) spinnerLabel_.getGraphics().dispose();
    spinnerLabel_ = null;
    if (cancelExport_!=null) cancelExport_.getGraphics().dispose();
    cancelExport_ = null;
    if (attentionProbes_!=null) attentionProbes_.clear();
    attentionProbes_ = null;
    if (selectedMolecules_!=null) selectedMolecules_.clear();
    selectedMolecules_ = null;
    if (applySettingsPopup_!=null) applySettingsPopup_.removeAll();
    applySettingsPopup_ = null;
    if (removeAnalytePopup_!=null) removeAnalytePopup_.removeAll();
    removeAnalytePopup_ = null;
    selectItem_ = null;
    deselectItem_ = null;
    deselectItem_ = null;
    applyToAllDoubles_ = null;
    lastClickedAnalyte_ = null;
    if (modifications_!=null) modifications_.clear();
    modifications_ = null;
    if (chromExport_!=null && chromExport_.getGraphics()!=null) chromExport_.getGraphics().dispose();
    chromExport_ = null;
    if (fromShortToExpName_!=null) fromShortToExpName_.clear();
    fromShortToExpName_ = null;
    rtTolerance_ = null;
    if (chromExportThread_!=null) chromExportThread_.cleanup();
    chromExportThread_ = null;
    if (timer_!=null){
      timer_.cancel();
      timer_.purge();
    }
    timer_ = null;
  }
  
  public Hashtable<String,Hashtable<String,Vector<Double>>> getResultValues() throws NumberFormatException, CalculationNotPossibleException{
    int exportType = ExportOptionsVO.EXPORT_NO_DEVIATION;
    if (isGrouped_) exportType = ExportOptionsVO.EXPORT_SD_DEV_AND_ERROR;
    ExportOptionsVO expOptions = new ExportOptionsVO(exportType,"1", false, true, false, 0);
    ResultDisplaySettingsVO settings = new ResultDisplaySettingsVO(settingsVO_);
    //TODO: this is set since I do not know about the output options of mzTab - maybe use the real settingsVO_ in future
    settings.setType("relative value");
    return HeatMapDrawing.extractValuesOfInterest(resultsOfOneGroup_, Integer.parseInt((String)maxIsotopes_.getSelectedItem()), settings, null, expOptions,modifications_);
  }
}
