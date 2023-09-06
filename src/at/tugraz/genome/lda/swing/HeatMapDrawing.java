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
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.Set;
import java.util.Timer;
import java.util.TimerTask;
import java.util.Vector;
import java.util.concurrent.ConcurrentHashMap;

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

import at.tugraz.genome.dbutilities.SimpleValueObject;
import at.tugraz.genome.lda.ChromExportThread;
import at.tugraz.genome.lda.LipidDataAnalyzer;
import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.TooltipTexts;
import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.analysis.ComparativeAnalysis;
import at.tugraz.genome.lda.analysis.HeatMapClickListener;
import at.tugraz.genome.lda.analysis.LipidomicsHeatMap;
import at.tugraz.genome.lda.analysis.exception.CalculationNotPossibleException;
import at.tugraz.genome.lda.exception.ExcelInputFileException;
import at.tugraz.genome.lda.exception.ExportException;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.exception.NoRuleException;
import at.tugraz.genome.lda.exception.RetentionTimeGroupingException;
import at.tugraz.genome.lda.exception.RulesException;
import at.tugraz.genome.lda.export.ExcelAndTextExporter;
import at.tugraz.genome.lda.msn.RulesContainer;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.AutoAnalyteAddVO;
import at.tugraz.genome.lda.vos.ExportOptionsVO;
import at.tugraz.genome.lda.vos.ResultAreaVO;
import at.tugraz.genome.lda.vos.ResultCompVO;
import at.tugraz.genome.lda.vos.ResultDisplaySettingsVO;
import at.tugraz.genome.maspectras.parser.exceptions.SpectrummillParserException;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class HeatMapDrawing extends JPanel implements ActionListener
{
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
//  protected OmegaExportDialog omegaExport_;
  protected ResultSelectionSettings combinedChartSettings_;
  protected boolean isGrouped_ = false;
  private JFileChooser exportFileChooser_;
  protected JLabel exportLabel_;
  protected JProgressBar exportProgress_;
  protected JLabel spinnerLabel_;
  protected JButton cancelExport_;
  
  private Set<Integer> selectedMoleculeRows_;
  /** when single items in the heat map are selected; first key: analyte name; second key: experiment name*/
  private Hashtable<String,Hashtable<String,Color>> selectedSingleMolecules_;
  /** stores the possible modifications when single items in the heat map are selected; first key: analyte name; second key: modification; value: modification*/
  private Hashtable<String,Hashtable<String,String>> selectedSingleMoleculesMods_;
  /** absolute file paths for the experiments where single items in the heat map are selected; first key: analyte name; value: absolute file path*/  
  private Hashtable<String,String> selectedSingleMoleculesAbsPaths_;
  /** stores information about the experiments where the automated quantitation should take place when single items in the heat map are selected;
   * first key: analyte name; value: vector of objects holding information about the other experiments*/
  private Hashtable<String,Vector<AutoAnalyteAddVO>> selectedSingleMoleculesAutoAnalyteAddVO_;
  /** stores the highest isotope number when single items in the heat map are selected; first key: analyte name; second key: experiment name*/
  private Hashtable<String,Integer> selectedSingleMoleculesMaxIso_;
  
  private boolean combinedDialogOpen_ = false;
  private boolean parentAction_;
  
  private PopupMenu applySettingsPopup_;
  private PopupMenu removeAnalytePopup_;
  
  private MenuItem selectItem_;
  private MenuItem deselectItem_;

  /** option "Select" in the popup menu appearing when a cell in the heat map is clicked by the right mouse button*/
  private MenuItem selectSingleItem_;
  /** option "Deselect" in the popup menu appearing when a cell in the heat map is clicked by the right mouse button*/
  private MenuItem deselectSingleItem_;

  
  private MenuItem applyToAllDoubles_;
//  private String lastClickedExp_;
//  private String lastClickedMol_;
//  private ResultCompVO lastClickedResultVO_;
  private int[] lastClickedCellPos_;
//  private String lastClickedAnalyte_;
  private int lastClickedRow_;
  private ArrayList<String> modifications_;
  private ChromExportDialog chromExport_;
  private Hashtable<String,String> fromShortToExpName_;
  private Double rtTolerance_;
  private HeatMapDrawing ungroupedPartner_;
  
  private ChromExportThread chromExportThread_;
  private Timer timer_;
  
  
  /**
   * constructor for a heat map
   * @param resultsOfOneGroup 					the values; first key molecule name; second key: experiment name
   * @param experimentNames 						the names of the experiments/sample groups
   * @param moleculeNames 							the names of the analytes
   * @param isLookup 										the lookup numbers for the individual internal standard correction options; key: standard correction option; value: lookup number
   * @param esLookup 										the lookup numbers for the individual external standard correction options; key: standard correction option; value: lookup number
   * @param statusText 									label containing information about the lipid the mouse is currently hoverd over (displayed at bottom of heat map)
   * @param listener 										call back listener
   * @param groupName 									analyte class name
   * @param ungroupedPartner 						for sample groups, the heat map containing the individual experiments, null otherwise
   * @param displaySettings 						dialog box for choosing the displayed values
   * @param selectionSettings 					a dialog box for choosing the values to be displayed in the heat map
   * @param combinedChartSettings 			a dialog box for choosing the values to be displayed in a combined chart
   * @param exportSettings 							a dialog box for choosing parameters for the export
   * @param analysisModule							the analysis module containing the params and data about max isotopes, modifications, rt tolerance
   */
  public HeatMapDrawing(Hashtable<String,Hashtable<String,ResultCompVO>> resultsOfOneGroup, Vector<String> experimentNames,
  		Vector<String> moleculeNames, Hashtable<String,Integer> isLookup, Hashtable<String,Integer> esLookup,
  		JLabel statusText, HeatMapClickListener listener, String groupName, HeatMapDrawing ungroupedPartner, 
      ResultDisplaySettings displaySettings, ResultSelectionSettings selectionSettings, 
      ResultSelectionSettings combinedChartSettings, ExportSettingsPanel exportSettings, ComparativeAnalysis analysisModule)
  {
  	this.resultsOfOneGroup_ = resultsOfOneGroup;
    this.experimentNames_ = new Vector<String>(experimentNames);
    this.moleculeNames_ = moleculeNames;
    this.isLookup_ = isLookup;
    this.esLookup_ = esLookup;
    this.statusText_ = statusText;
    this.heatMapListener_ = listener;
    this.groupName_ = groupName;
    this.isGrouped_ = ungroupedPartner!=null ? true : false;
    this.ungroupedPartner_ = ungroupedPartner;
    this.displaySettings_ = displaySettings;
    this.displaySettings_.addActionListener(this);
    this.selectionSettings_ = selectionSettings;
    this.selectionSettings_.addActionListener(this);
    this.combinedChartSettings_ = combinedChartSettings;
    this.combinedChartSettings_.addActionListener(this);
    this.exportSettings_ = exportSettings;
    this.maxIsotopes_ = generateMaxIsotopesJComboBox(analysisModule.getMaxIsotopesOfGroup(groupName));
    this.modifications_ = new ArrayList<String>(analysisModule.getModifications().get(groupName));
    this.rtTolerance_ = analysisModule.getRtTolerance();
    
    
    settingsVO_ = displaySettings_.getSettingsVO();
    
    selectedMoleculeRows_ = ConcurrentHashMap.newKeySet();
    selectedSingleMolecules_ = new Hashtable<String,Hashtable<String,Color>>();
    selectedSingleMoleculesMods_ = new Hashtable<String,Hashtable<String,String>>(); 
    selectedSingleMoleculesAbsPaths_ = new Hashtable<String,String>();
    selectedSingleMoleculesAutoAnalyteAddVO_ = new Hashtable<String,Vector<AutoAnalyteAddVO>>();
    selectedSingleMoleculesMaxIso_ = new Hashtable<String,Integer>();
    parentAction_ = true;
    
    combinedDialogOpen_ = false;
    
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
    selectSingleItem_ = new MenuItem("Select");
    selectSingleItem_.addActionListener(this);
    applySettingsPopup_.add(selectSingleItem_);
    deselectSingleItem_ = new MenuItem("Deselect");
    deselectSingleItem_.addActionListener(this);
    applySettingsPopup_.add(deselectSingleItem_);

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
    
    exportFileChooser_ = new JFileChooser();
    exportFileChooser_.setPreferredSize(new Dimension(600,500));
    
    this.generateHeatMap();
    
    if (!isGrouped_){     
      chromExport_ = new ChromExportDialog("Chrom export",getDisplayNames(),this.heatmap_.getAnalyteNames(),modifications_,this);
////      chromExport_.addActionListener(this);
    }
    chromExportThread_ = null;
    this.initTimer();
    
  }
  
  private JComboBox<String> generateMaxIsotopesJComboBox(int maxIsotopesOfGroup)
  {
  	JComboBox<String> jComboBox = new JComboBox<String>();
    for (int i=0;i!=maxIsotopesOfGroup;i++){
    	jComboBox.addItem(String.valueOf(i));
    }
    jComboBox.setSelectedItem(String.valueOf(maxIsotopesOfGroup-1));
    SelectionItemListener isotopeListener = new SelectionItemListener("ChangeIsotope");
    jComboBox.addItemListener(isotopeListener);
    jComboBox.setToolTipText(TooltipTexts.HEATMAP_ISOTOPES);
    return jComboBox;
  }
  
  
  public void generateHeatMap(){
    if (this.heatMapPanel_!=null)
      this.remove(heatMapPanel_);
    
    Vector<String> molNames = this.getSelectedMoleculeNames();
    
    try
    {
      int maxIsotopes = 0;
      try
      {
        maxIsotopes = Integer.parseInt((String)maxIsotopes_.getSelectedItem());
      } catch (NumberFormatException nfx){}  
      
      this.heatmap_ = new LipidomicsHeatMap(resultsOfOneGroup_, experimentNames_, heatMapListener_, molNames, maxIsotopes, settingsVO_, isMarkDoublePeaks());
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
    if (!isGrouped_ && chromExport_ != null) 
    	chromExport_.refreshNames(getDisplayNames());
  }
  
  
  public Dimension getTotalSize() {
    return new Dimension(renderedImage_.getWidth()+10,renderedImage_.getHeight());
  }
  
  public void actionPerformed(ActionEvent event){
    actionPerformed(event.getActionCommand());
  }
  
  @SuppressWarnings("unchecked")
  private void actionPerformed(String actionCommand)
  {
    Vector<String> foundUpdateables = new Vector<String>();
    Vector<AutoAnalyteAddVO> updateableAndAnalyteBefore = new Vector<AutoAnalyteAddVO>();
    Vector<SimpleValueObject> templateSpecies = new Vector<SimpleValueObject>();
    Hashtable<String,String> modHash = new Hashtable<String,String>();
    ResultCompVO vo = null;
    int maxIsotope = 0;
    if (actionCommand.equalsIgnoreCase(CHANGE_IS_STATUS)||
        actionCommand.equalsIgnoreCase(CHANGE_ES_STATUS)||
        actionCommand.equalsIgnoreCase(CHANGE_DOUBLE_STATUS)||
        actionCommand.equalsIgnoreCase("AcceptDisplaySettings")){
      settingsVO_ = displaySettings_.getSettingsVO();
      selectionSettings_.setVisible(false);
      boolean update = true;    
      if (combinedDialogOpen_){
        update = false;
        if (combinedChartSettings_.getSelected().size()==0){
          new WarningMessage(new JFrame(), "Error", "You have to select at least one molecule!");
          return;
        }else{
          combinedChartSettings_.setVisible(false);
          combinedDialogOpen_ = false;
          int maxIsotopes = Integer.parseInt((String)maxIsotopes_.getSelectedItem());
          Hashtable<String,String> preferredUnits = new Hashtable<String,String>();
          for (String name : combinedChartSettings_.getSelected())
          {
          	int index = heatmap_.getAnalyteIndex(name);
          	if (index >= 0)
          	{
          		preferredUnits.put(name, heatmap_.getPreferredUnit(index));
          	}
          }
          String preferredUnit = heatmap_.extractPreferredUnitForExp(preferredUnits);
          if (isGrouped_){
            heatMapListener_.combinedAnalyteGroupSelected(combinedChartSettings_.getSelected(), groupName_, maxIsotopes, rtTolerance_!=null,
            settingsVO_, preferredUnit, StaticUtils.getCorrespondingUnit(settingsVO_, preferredUnit,true));
          }else
            heatMapListener_.combinedAnalyteSelected(combinedChartSettings_.getSelected(), groupName_, maxIsotopes, rtTolerance_!=null,
            settingsVO_, preferredUnit, StaticUtils.getCorrespondingUnit(settingsVO_, preferredUnit,true));
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
        heatMapListener_.changeISStatus(groupName_, isGrouped_, showInternalStandards_.isSelected());
      if (actionCommand.equalsIgnoreCase(CHANGE_ES_STATUS) && parentAction_)
        heatMapListener_.changeESStatus(groupName_, isGrouped_, showExternalStandards_.isSelected());
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
    } else if (actionCommand.equalsIgnoreCase("openCombinedDialog")){
      combinedChartSettings_.setVisible(true);
      combinedDialogOpen_ = true;
    } else if (actionCommand.equalsIgnoreCase(ExportPanel.EXPORT_PNG)){
      exportFileChooser_.setFileFilter(new FileNameExtensionFilter("PNG (*.png)","png"));
      int returnVal = exportFileChooser_.showSaveDialog(new JFrame());
      if (returnVal == JFileChooser.APPROVE_OPTION) {
        File fileToStore = exportFileChooser_.getSelectedFile();
        @SuppressWarnings("rawtypes")
        Vector results = checkFileStorage(fileToStore,"png",this);
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
        Vector results = checkFileStorage(fileToStore,"svg",this);
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
      ExportOptionsVO expOptions = getExportOptions();
      LinkedHashMap<String,String> expFullPaths = heatMapListener_.getSampleResultFullPaths();
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
            Vector results = checkFileStorage(fileToStore,"xlsx",this);
            fileToStore = (File)results.get(0);
            if ((Boolean)results.get(1)){
              try {
                maxIsotope = Integer.parseInt((String)maxIsotopes_.getSelectedItem());
                Hashtable<String,Hashtable<String,ResultCompVO>> compVOs = resultsOfOneGroup_;
                if (isGrouped_)
                  compVOs = ungroupedPartner_.resultsOfOneGroup_;
                String preferredUnit = heatmap_.extractPreferredUnitForExp();
                Hashtable<String,Hashtable<String,Vector<Double>>> resultValues = HeatMapDrawing.extractValuesOfInterest(compVOs, maxIsotope, settingsVO_, preferredUnit, expOptions,modifications_);
                preferredUnit = StaticUtils.getCorrespondingUnit(settingsVO_,preferredUnit,true);
                boolean exportDoubleBondPositionsForClass = expOptions.isExportDoubleBondPositions();
                short speciesType = expOptions.getSpeciesType();
                if (exportDoubleBondPositionsForClass) {
                  int numberOfChains = 2;
                  try {
                    numberOfChains = Integer.parseInt(RulesContainer.getAmountOfChains(StaticUtils.getRuleName(groupName_, modifications_.get(0))));
                  } catch (RulesException | NoRuleException | IOException | SpectrummillParserException ex) {
                    System.out.println(ex.getMessage());
                  }
                  if (numberOfChains > 1 && speciesType == LipidomicsConstants.EXPORT_ANALYTE_TYPE_SPECIES) {
                    exportDoubleBondPositionsForClass = false;
                  //for lipid species with only one FA chain we export the double bond position information on species level
                  } else if (numberOfChains == 1 && speciesType != LipidomicsConstants.EXPORT_ANALYTE_TYPE_SPECIES) {
                    speciesType = LipidomicsConstants.EXPORT_ANALYTE_TYPE_SPECIES;
                  }
                }
                ExcelAndTextExporter.exportToFile(true,speciesType, exportDoubleBondPositionsForClass, groupName_,new BufferedOutputStream(new FileOutputStream(fileToStore)),true,
                    maxIsotope, getSelectedMoleculeNames(), rtTolerance_!=null, isGrouped_, experimentNames_,expIdToString, expFullPaths, heatMapListener_.getSamplesOfGroups(),
                    resultValues, preferredUnit, StaticUtils.getAreaTypeString(settingsVO_), expOptions, heatMapListener_.getComparativeResultsLookup(), modifications_);
              } 
              catch (NumberFormatException e) {new WarningMessage(new JFrame(), "Error", e.getMessage());}
              catch (FileNotFoundException e) {new WarningMessage(new JFrame(), "Error", e.getMessage());}
              catch (IOException e) {new WarningMessage(new JFrame(), "Error", e.getMessage());}
              catch (CalculationNotPossibleException | ExcelInputFileException | ExportException | SpectrummillParserException | LipidCombinameEncodingException | RetentionTimeGroupingException e) {
                new WarningMessage(new JFrame(), "Error", e.getMessage());
              }
            }  
          }
        } else if (actionCommand.equalsIgnoreCase(ExportPanel.EXPORT_TEXT)){
          exportFileChooser_.setFileFilter(new FileNameExtensionFilter("Text (*.txt)","txt"));
          int returnVal = exportFileChooser_.showSaveDialog(new JFrame());
          if (returnVal == JFileChooser.APPROVE_OPTION) {
            File fileToStore = exportFileChooser_.getSelectedFile();
            @SuppressWarnings("rawtypes")
            Vector results = checkFileStorage(fileToStore,"txt",this);
            fileToStore = (File)results.get(0);
            if ((Boolean)results.get(1)){
              try {
                maxIsotope = Integer.parseInt((String)maxIsotopes_.getSelectedItem());
                Hashtable<String,Hashtable<String,ResultCompVO>> compVOs = resultsOfOneGroup_;
                if (isGrouped_)
                  compVOs = ungroupedPartner_.resultsOfOneGroup_;
                String preferredUnit = heatmap_.extractPreferredUnitForExp();
                Hashtable<String,Hashtable<String,Vector<Double>>> resultValues = HeatMapDrawing.extractValuesOfInterest(compVOs, maxIsotope, settingsVO_, preferredUnit, expOptions,modifications_);
                preferredUnit = StaticUtils.getCorrespondingUnit(settingsVO_,preferredUnit,true);
                boolean exportDoubleBondPositionsForClass = expOptions.isExportDoubleBondPositions();
                short speciesType = expOptions.getSpeciesType();
                if (exportDoubleBondPositionsForClass) {
                  int numberOfChains = 2;
                  try {
                    numberOfChains = Integer.parseInt(RulesContainer.getAmountOfChains(StaticUtils.getRuleName(groupName_, modifications_.get(0))));
                  } catch (RulesException | NoRuleException | IOException | SpectrummillParserException ex) {
                    System.out.println(ex.getMessage());
                  }
                  if (numberOfChains > 1 && speciesType == LipidomicsConstants.EXPORT_ANALYTE_TYPE_SPECIES) {
                    exportDoubleBondPositionsForClass = false;
                  //for lipid species with only one FA chain we export the double bond position information on species level
                  } else if (numberOfChains == 1 && speciesType != LipidomicsConstants.EXPORT_ANALYTE_TYPE_SPECIES) {
                    speciesType = LipidomicsConstants.EXPORT_ANALYTE_TYPE_SPECIES;
                  }
                }
                ExcelAndTextExporter.exportToFile(true,speciesType, exportDoubleBondPositionsForClass, groupName_,new BufferedOutputStream(new FileOutputStream(fileToStore)),false,
                    maxIsotope, getSelectedMoleculeNames(), rtTolerance_!=null, isGrouped_, experimentNames_,expIdToString, expFullPaths, heatMapListener_.getSamplesOfGroups(),
                    resultValues, preferredUnit, StaticUtils.getAreaTypeString(settingsVO_), expOptions, heatMapListener_.getComparativeResultsLookup(),
                    modifications_);

              }
              catch (NumberFormatException e) {new WarningMessage(new JFrame(), "Error", e.getMessage());}
              catch (FileNotFoundException e) {new WarningMessage(new JFrame(), "Error", e.getMessage());}
              catch (IOException e) {new WarningMessage(new JFrame(), "Error", e.getMessage());}
              catch (CalculationNotPossibleException | ExcelInputFileException | ExportException | SpectrummillParserException | LipidCombinameEncodingException | RetentionTimeGroupingException e) {
                new WarningMessage(new JFrame(), "Error", e.getMessage());
              }
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
    	    ExportOptionsVO expOptions = getExportOptions();
    	    File exportFile = exportFileChooser_.getSelectedFile();
    	    if (actionCommand.equalsIgnoreCase(ExportPanel.EXPORT_MZTAB)){
          File fileToStore = exportFileChooser_.getSelectedFile();
          @SuppressWarnings("rawtypes")
          Vector results = checkFileStorage(fileToStore,"txt",this);
          fileToStore = (File)results.get(0);
          if ((Boolean)results.get(1))
    	        heatMapListener_.exportMzTab(fileToStore, expOptions.getSpeciesType(), expOptions.isExportDoubleBondPositions());
    	    }
    	    else if (actionCommand.equalsIgnoreCase(ExportPanel.EXPORT_RDB)) heatMapListener_.exportRdb(exportFile);
      }
      exportFileChooser_.setSelectedFile(new File(""));
    	
    } else if (actionCommand.equalsIgnoreCase(ExportPanel.EXPORT_CHROMS)){
      chromExport_.setVisible(true);
    } 
    else if (actionCommand.equalsIgnoreCase(ExportPanel.EXPORT_MAF)) {
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
      boolean chromFileExists = false;
      for (int i=0;i!=this.experimentNames_.size();i++ ){
        String name = this.experimentNames_.get(i);
        if (chromExport_.isExperimentSelected(heatMapListener_.getDisplayName(name))){
        	String resultFile = heatmap_.getAbsoluteFilePathOfExperiment(i);
        	if (resultFile != null)
        	{
        		resultFiles.add(resultFile);
            String chromFileBase = StaticUtils.extractChromBaseName(resultFile,name);
            if (chromFileBase!=null && chromFileBase.length()>0){
              chromFileExists = true;        
            }
            if (chromFileExists && !StaticUtils.existChromFiles(chromFileBase)) chromFileExists = false;           
            if (chromFileExists){
              String chromFilePath = chromFileBase+".chrom";
              chromsToUse.add(chromFilePath);
              expsToExport.add(heatMapListener_.getDisplayName(name));
            } else {
              new WarningMessage(new JFrame(), "Error", "The '.chrom' files are required to be in the same directory as the specified result excel files!");
              break;
            }
        	} 
        }
      } if (chromFileExists) {
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
              results = checkFileStorage(fileToStore,"png",this);
            if (chromExport_.getPictureType().equalsIgnoreCase(ExportPanel.EXPORT_SVG))
              results = checkFileStorage(fileToStore,"svg",this);
            fileToStore = (File)results.get(0);
            if ((Boolean)results.get(1) && chromExportThread_==null){
              chromExportThread_ = new ChromExportThread(groupName_,fileToStore,chromExport_.getPictureDimension(),
                  expsToExport,resultFiles,chromsToUse,analsToExport,chromExport_.getPictureType(),exportSettings_.getSettings().isAnalyteInColumn(),
                  rtTolerance_,isLookup_,esLookup_);
              chromExportThread_.start();
              spinnerLabel_.setVisible(true);
              exportProgressPanel_.setVisible(true);
            } else if (chromExportThread_!=null)
            {
            	new WarningMessage(new JFrame(), "Error", "There is already a chromatogram export running! Wait until it is finished!");
            	return;
            }
          }
        }
      }
    } else if (actionCommand.equalsIgnoreCase("Choose just one peak for doubles") || actionCommand.equalsIgnoreCase("Quant. anal. at not found")||
        actionCommand.equalsIgnoreCase("Take exact peak for others") || actionCommand.equalsIgnoreCase("Select"))
    {
      String expName = experimentNames_.get(lastClickedCellPos_[0]);
      String molName = heatmap_.getOriginalAnalyteName(lastClickedCellPos_[1]);
      if ((actionCommand.equalsIgnoreCase("Quant. anal. at not found")||actionCommand.equalsIgnoreCase("Take exact peak for others")) &&
          selectedSingleMolecules_.containsKey(molName) && !selectedSingleMolecules_.get(molName).keySet().iterator().next().equalsIgnoreCase(expName))
      {
        new WarningMessage(new JFrame(), "Error", "It is not allowed to select more than one from the same lipid species");
        return;
      }
      vo = heatmap_.getCompVO(lastClickedCellPos_[0],lastClickedCellPos_[1]);
      int sumCompRow = heatmap_.getSumCompVORow(lastClickedCellPos_[1]);
      Hashtable<String,String> availableMods = new Hashtable<String,String>();
      for (String modName : modifications_){
        if (vo.containsMod(modName)) availableMods.put(modName, modName);
      }
      int maxIsotopes = Integer.parseInt((String)maxIsotopes_.getSelectedItem());
      foundUpdateables = new Vector<String>();
      updateableAndAnalyteBefore = new Vector<AutoAnalyteAddVO>();
      modHash = new Hashtable<String,String>();
      maxIsotope = 0;
      for (int i=0;i!=experimentNames_.size();i++){
        ResultCompVO otherVO = heatmap_.getCompVO(i, sumCompRow);
        if (actionCommand.equalsIgnoreCase("Choose just one peak for doubles")){
          if (i!=lastClickedCellPos_[0] && otherVO.getMoreThanOnePeak(otherVO.getAvailableIsotopeNr(maxIsotopes))){
            boolean foundMod = false;
            for (String modName : availableMods.keySet()){
              if (otherVO.getMoreThanOnePeak(otherVO.getAvailableIsotopeNr(maxIsotopes), modName)){
                modHash.put(modName, modName);
                foundMod = true;
              }
            }
            if (foundMod)
              foundUpdateables.add(otherVO.getAbsoluteFilePath());
          }
        } else if (actionCommand.equalsIgnoreCase("Quant. anal. at not found")||actionCommand.equalsIgnoreCase("Take exact peak for others")||
            actionCommand.equalsIgnoreCase("Select")){
          if (i==lastClickedCellPos_[0])
            continue;
          boolean modMissing = false;
          if (!otherVO.existsInFile())
            modMissing = true;
          else{
            for (String modName : availableMods.keySet()){
              if (!otherVO.containsMod(modName)){
                modMissing = true;
                break;
              }
            }
          }
          if (modMissing){
            String previousElement = "";
            for (int j=sumCompRow-1;j>-1;j--){
              if (heatmap_.getSumCompVO(i,j).existsInFile()){
                previousElement = heatmap_.getOriginalAnalyteName(j);
                break;
              }
            }
            updateableAndAnalyteBefore.add(new AutoAnalyteAddVO(experimentNames_.get(i),otherVO.getAbsoluteFilePath(), previousElement));
            maxIsotope = vo.getUsedIsotpes();
            for (String modName : modifications_){
              if (vo.containsMod(modName))modHash.put(modName, modName);
            }
          }
        }
      }
      if (foundUpdateables.size()>0 || updateableAndAnalyteBefore.size()>0){
        templateSpecies.add(new SimpleValueObject(molName,expName));
      }else{  
        if (actionCommand.equalsIgnoreCase("Choose just one peak for doubles"))
        {
        	new WarningMessage(new JFrame(), "Error", "There are no double peaks to eliminate for "+molName+"!");
        	return;
        }
        else if (actionCommand.equalsIgnoreCase("Quant. anal. at not found") || actionCommand.equalsIgnoreCase("Take exact peak for others") ||
            actionCommand.equalsIgnoreCase("Select"))
        {
        	new WarningMessage(new JFrame(), "Error", "There are no peaks to add for "+molName+"!");
        	return;
        }
      }
    } else if (actionCommand.equalsIgnoreCase("Remove analyte in all probes")){
    	Vector<String> messages = new Vector<String>();
    	Set<Integer> analytesToRemove = ConcurrentHashMap.newKeySet();
    	analytesToRemove.addAll(selectedMoleculeRows_);
    	analytesToRemove.add(lastClickedRow_);
    	generateAnalyteRemovalDialogue(analytesToRemove, messages);
    	
    	Set<Integer> analytesToModify = ConcurrentHashMap.newKeySet();
    	for (int row : analytesToRemove)
    	{
    		if (heatmap_.isMolecularSpeciesLevel(row))
    		{
    			analytesToModify.add(row);
    		}
    	}
    	analytesToRemove.removeAll(analytesToModify);
    	
      Vector<String> selectedMods = CheckBoxOptionPane.showConfirmDialog(new JFrame(), "Confirmation", messages, modifications_,false);
      Set<String> filePaths = heatmap_.eliminateMolecularSpecies(analytesToModify, selectedMods);
      
      if (selectedMods.size()>0)
      {
      	Set<ResultCompVO> toRemove = heatmap_.getPresentCompVOsWithMod(analytesToRemove, selectedMods);
        heatMapListener_.eliminateAnalyteEverywhere(groupName_, toRemove, selectedMods, filePaths);
      }
    }else if (actionCommand.equalsIgnoreCase("Select analyte") || actionCommand.equalsIgnoreCase("Deselect analyte")){  
      Graphics2D g2 = (Graphics2D)renderedImage_.getGraphics();
      this.selectAnalyte(lastClickedRow_);
      if (actionCommand.equalsIgnoreCase("Select analyte")){
        this.selectAnalyte(lastClickedRow_);
        g2.setColor(Color.BLUE);
      }else if (actionCommand.equalsIgnoreCase("Deselect analyte")){
        g2.setColor(Color.BLACK);
        this.deselectAnalyte(lastClickedRow_);
      }  
      Rectangle rectForName = heatmap_.getRectangleForRowNumber(lastClickedRow_);
      rectForName.setLocation(new Point(rectForName.getLocation().x+imagePositionX_,rectForName.getLocation().y+imagePositionY_));
      g2.fillRect(rectForName.x+1,rectForName.y+1, rectForName.width-2, rectForName.height-2);
      Font descriptionFont = new Font("Dialog",Font.PLAIN, 9);
      FontMetrics descriptionFontMetrics = g2.getFontMetrics();
      int textHeight = descriptionFontMetrics.getHeight();
      g2.setColor(Color.WHITE);
      g2.setFont(descriptionFont);
      g2.drawString(heatmap_.getAnalyteName(lastClickedRow_),rectForName.x+1,rectForName.y + textHeight/2);
      this.invalidate();
      this.updateUI();
    } else if (actionCommand.equalsIgnoreCase("stopChromExport")){
      chromExportThread_.stopThread();
    } else if (actionCommand.equalsIgnoreCase("exportSelectionDialog")){
      heatMapListener_.showExportSettingsDialog(this.isGrouped_);
    }
    if ((foundUpdateables.size()>0 || updateableAndAnalyteBefore.size()>0) && actionCommand.equalsIgnoreCase("Choose just one peak for doubles") || actionCommand.equalsIgnoreCase("Quant. anal. at not found")||
        actionCommand.equalsIgnoreCase("Take exact peak for others")){
      SimpleValueObject currentlySelected = templateSpecies.get(0);
      @SuppressWarnings("rawtypes")
      Vector result = getSelectedIndividualSpeciesInCorrectOrder(currentlySelected,modHash,vo.getAbsoluteFilePath(),updateableAndAnalyteBefore,maxIsotope);
      templateSpecies = (Vector<SimpleValueObject>)result.get(0);
      Hashtable<String,String> allMods = (Hashtable<String,String>)result.get(1);
      Vector<String> absPaths = (Vector<String>)result.get(2);
      Hashtable<String,Vector<AutoAnalyteAddVO>> updateablesAndAnalytesBefore = (Hashtable<String,Vector<AutoAnalyteAddVO>>)result.get(3);
      Vector<Integer> maxIsotopes = (Vector<Integer>)result.get(4);
      
      Vector<String> confirmationMessages = new Vector<String>();
      String confirmationMessage = "Are you sure you want to eliminate double peaks for "+currentlySelected.getLabel()+"? The peak nearest to the RT of "+currentlySelected.getValue()+" is taken!";
      if (actionCommand.equalsIgnoreCase("Quant. anal. at not found")||actionCommand.equalsIgnoreCase("Take exact peak for others")){
        confirmationMessage = "Are you sure you want to automatically add peaks for the following molecules? ";
        if (actionCommand.equalsIgnoreCase("Quant. anal. at not found"))
          confirmationMessage += "The peak nearest to the RT of the following experiments are taken:";
        else if (actionCommand.equalsIgnoreCase("Take exact peak for others"))
          confirmationMessage += "The exact peak borders of the following experiments are taekn:";
        confirmationMessages.add(confirmationMessage);
        for (SimpleValueObject species : templateSpecies)
          confirmationMessages.add(species.getLabel()+":    "+species.getValue());
      }else{
        confirmationMessages.add(confirmationMessage);
      }
      Vector<String> selectedMods = CheckBoxOptionPane.showConfirmDialog(new JFrame(), "Confirmation", confirmationMessages, new ArrayList<String>(allMods.values()),true);
      if (selectedMods.size()>0){
        boolean exactProbePosition = false;
        if (actionCommand.equalsIgnoreCase("Take exact peak for others")) exactProbePosition = true;
        if (actionCommand.equalsIgnoreCase("Choose just one peak for doubles"))
          heatMapListener_.eliminateDoublePeaks(groupName_,templateSpecies.get(0).getLabel(),vo.getAbsoluteFilePath(),selectedMods,foundUpdateables);
        else if (actionCommand.equalsIgnoreCase("Quant. anal. at not found")||actionCommand.equalsIgnoreCase("Take exact peak for others")){
          Vector<String> molNames = new Vector<String>();
          for (SimpleValueObject svo : templateSpecies) molNames.add(svo.getLabel());
          heatMapListener_.addAnalytesEverywhereAtPosition(groupName_,molNames,absPaths,selectedMods,updateablesAndAnalytesBefore,maxIsotopes,exactProbePosition);          
        }
      }

    
    } else if ((actionCommand.equalsIgnoreCase("Select") && ((foundUpdateables.size()>0 || updateableAndAnalyteBefore.size()>0))) || actionCommand.equalsIgnoreCase("Deselect")){  
      Graphics2D g2 = (Graphics2D)renderedImage_.getGraphics();
      boolean printAttentionRectangle = false;
      if (actionCommand.equalsIgnoreCase("Select")){
        if (!this.selectAnalyte(moleculeNames_.get(lastClickedCellPos_[1]),experimentNames_.get(lastClickedCellPos_[0]),heatmap_.getColorForCell(lastClickedCellPos_[1], lastClickedCellPos_[0]),modHash,
            vo.getAbsoluteFilePath(),updateableAndAnalyteBefore,maxIsotope))
          return;
        g2.setColor(Color.BLUE);
      }else if (actionCommand.equalsIgnoreCase("Deselect")){
        g2.setColor(this.deselectAnalyte(moleculeNames_.get(lastClickedCellPos_[1]),experimentNames_.get(lastClickedCellPos_[0])));
        if (isMarkDoublePeaks() && heatmap_.getAttentionProbe(lastClickedCellPos_[1], experimentNames_.get(lastClickedCellPos_[0])) != null)
        {
          printAttentionRectangle = true;
        }
      }
      Rectangle cellRect = heatmap_.getRectangleForCellByRowAndColumn(lastClickedCellPos_[1], lastClickedCellPos_[0]);
      cellRect.setLocation(new Point(cellRect.getLocation().x+imagePositionX_,cellRect.getLocation().y+imagePositionY_));
      g2.fillRect(cellRect.x+1,cellRect.y+1,cellRect.width-1,cellRect.height-1);
      if (printAttentionRectangle){
        heatmap_.paintAttentionRectangle(g2, heatmap_.getAttentionProbe(moleculeNames_.get(lastClickedCellPos_[1]), experimentNames_.get(lastClickedCellPos_[0])), 
            lastClickedCellPos_[1], lastClickedCellPos_[0], heatmap_.getExpressionImageXStart(), heatmap_.getExpressionImageYStart());
      }
      this.invalidate();
      this.updateUI();
    }

  }
  
  private void generateAnalyteRemovalDialogue(Set<Integer> analytesToRemove, Vector<String> messages)
  {
    String analyteList = "";
    for (Integer analyteRow : analytesToRemove)
    {
    	if (heatmap_.isMolecularSpeciesLevel(analyteRow))
    	{
    		analyteList += String.format("%s (RT=%s min) <br>", heatmap_.getMolecularSpeciesLevelName(analyteRow), heatmap_.getRetentionTime(analyteRow));
    	}
    	else
    	{
    		ArrayList<Integer> rows = heatmap_.getAllRowsOfSameSumComposition(analyteRow);
    		analyteList += String.format("%s (RT=%s min), this sum composition contains %s identifications at the molecular species level, which will be deleted as well", 
    				heatmap_.getMolecularSpeciesLevelName(analyteRow), heatmap_.getRetentionTime(analyteRow), rows.size());
    		if (rows.size() > 0)
    		{
    			analyteList += ": ";
    		}
    		else
    		{
    			analyteList += ". <br>";
    		}
    		for (int i=0; i<rows.size(); i++)
    		{
    			int row = rows.get(i);
    			analyteList += String.format("%s", heatmap_.getMolecularSpeciesLevelName(row), heatmap_.getRetentionTime(row));
    			if (i < rows.size()-1)
    			{
    				analyteList += ", ";
    			}
    			else
    			{
    				analyteList += ". <br>";
    			}
    		}
    	}
    }
    messages.add(String.format("<html>Do you really want to delete the following identifications of the lipid class '%s' in all probes: <br>%s </html>", 
    		groupName_, analyteList));
  }
  
  private static Hashtable<String,Hashtable<String,Vector<Double>>> extractValuesOfInterest(Hashtable<String,Hashtable<String,ResultCompVO>> vos, int maxIsotope, ResultDisplaySettingsVO settingVO, String preferredUnit, ExportOptionsVO expOptions, ArrayList<String> modifications) throws CalculationNotPossibleException{
    Hashtable<String,Hashtable<String,Vector<Double>>> results = new Hashtable<String,Hashtable<String,Vector<Double>>>();
    for (String molKey : vos.keySet()){
      Hashtable<String,Vector<Double>> expValues = new Hashtable<String,Vector<Double>>();
      for (String expKey : vos.get(molKey).keySet()){
        ResultCompVO compVO = vos.get(molKey).get(expKey);
        int isoNr = maxIsotope;
        if (!settingVO.getType().equalsIgnoreCase(ResultDisplaySettingsVO.REL_MEASURED_CLASS_AMOUNT) && 
            !settingVO.getType().equalsIgnoreCase(ResultDisplaySettingsVO.REL_TOTAL_AMOUNT))
          isoNr = compVO.getAvailableIsotopeNr(maxIsotope);
        double myArea = compVO.getArea(isoNr, settingVO);
        if (settingVO.getType().equalsIgnoreCase(ResultDisplaySettingsVO.REL_MEASURED_CLASS_AMOUNT) || 
            settingVO.getType().equalsIgnoreCase(ResultDisplaySettingsVO.REL_TOTAL_AMOUNT))
          compVO.getArea(maxIsotope, settingVO);
        myArea = StaticUtils.getAreaInCorrespondingUnit(myArea, preferredUnit);
        Vector<Double> areaPlusDev = new Vector<Double>();
        areaPlusDev.add(myArea);
        expValues.put(expKey, areaPlusDev);
      }   
      results.put(molKey, expValues);
    }
    return results;
  }
  
  @SuppressWarnings({ "unchecked", "rawtypes" })
  public static Vector checkFileStorage (File file, String suffix, JPanel panel){
    Vector results = new Vector();
    boolean store = true;
    File fileToStore = new File(file.getAbsolutePath());
    if (fileToStore.exists()){
      if (JOptionPane.showConfirmDialog(panel, "The file "+fileToStore.getName()+" exists! Replace existing file?") != JOptionPane.YES_OPTION)
        store = false;
    }else{
      if (fileToStore.getName().indexOf(".")==-1){
        fileToStore = new File(fileToStore.getAbsoluteFile()+"."+suffix);
        if (fileToStore.exists()){
          if (JOptionPane.showConfirmDialog(panel, "The file "+fileToStore.getName()+" exists! Replace existing file?") != JOptionPane.YES_OPTION)
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
              heatmap_.getExpressionImageYStart()<=yInImage&&yInImage<heatmap_.getExpressionImageYEnd())
          {
          	int[] cellPos = heatmap_.getCellPosition(xInImage, yInImage);
            String experiment = heatmap_.getExperimentName(cellPos[0]);
            String analyte = heatmap_.getAnalyteName(cellPos[1]);
            if (experiment!=null&&analyte!=null)
            {
              if (statusText_!=null){
                ResultCompVO compVO = heatmap_.getCompVO(cellPos[0],cellPos[1]);
                String molName = heatmap_.getMolecularSpeciesName(cellPos[1]);
                Double molContribution = heatmap_.getMolecularSpeciesContributionOfAllMods(experiment, cellPos[1]);
                
                int maxIsotope = Integer.parseInt((String)maxIsotopes_.getSelectedItem());
                String statusText = "";
                try {
                  double relativeValue = compVO.getRelativeValue(compVO.getAvailableIsotopeNr(maxIsotope),settingsVO_,molName,heatmap_.getMolecularSpeciesContributionOfAllMods(experiment, cellPos[1]));
                  String sampleTypeText = "Sample";
                  if (isGrouped_) sampleTypeText = "Group";
                  statusText = "Lipid: "+analyte+"; "+sampleTypeText+": "+experiment+"; Relative: "+StaticUtils.extractDisplayValue(relativeValue);
                  double value = compVO.getArea(compVO.getAvailableIsotopeNr(maxIsotope), settingsVO_);
                  if (settingsVO_.getType().equalsIgnoreCase(ResultDisplaySettingsVO.REL_MEASURED_CLASS_AMOUNT) || 
                      settingsVO_.getType().equalsIgnoreCase(ResultDisplaySettingsVO.REL_TOTAL_AMOUNT))
                    value = compVO.getArea(maxIsotope, settingsVO_);
                  String selectedAreaString = StaticUtils.getAreaString(value*molContribution, settingsVO_, heatmap_.getPreferredUnit(cellPos[1]));
                  if (selectedAreaString!=null&&selectedAreaString.length()>0)
                    statusText += "; "+selectedAreaString;
                  double standValue = compVO.getStandardizedArea(compVO.getAvailableIsotopeNr(maxIsotope), settingsVO_.getISStandMethod(),  settingsVO_.getESStandMethod(), settingsVO_.considerDilution());
                  statusText += "; Stand-Area: "+standValue*molContribution;
                  statusText += "; Quant-Area: "+compVO.getOriginalArea(compVO.getAvailableIsotopeNr(maxIsotope))*molContribution+"; Isos: "+(compVO.getAvailableIsotopeNr(maxIsotope)+1);
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
          	int rowNumber = heatmap_.getRowNumber(x, y);
            rectToDraw_ = heatmap_.getRectangleForRowNumber(rowNumber);
            this.drawARectangle(g2);
          } else if (heatmap_.getExpressionImageXStart()<=xInImage&&xInImage<heatmap_.getExpressionImageXEnd()&&
              heatmap_.getColumnNameStart()<=yInImage&&yInImage<heatmap_.getColumnNameEnd())
          {
            String expName = heatmap_.getColumnName(x,y);
            if (expName!=null&&expName.length()>0)
            {
              rectToDraw_ = heatmap_.getRectangleForColumnName(expName);
              this.drawARectangle(g2);              
            }
          } 
//          else if (heatmap_.getLipidDoubleBondRowNameStart()<xInImage&& xInImage<heatmap_.getLipidDoubleBondRowNameEnd() &&
//              heatmap_.getExpressionImageYStart()<=yInImage&&yInImage<heatmap_.getExpressionImageYEnd()){
//            String[] doubleBondInfo = heatmap_.getLipidDoubleBondRowName(x,y);
//            if (doubleBondInfo[1] != "") {
//              rectToDraw_ = heatmap_.getRectangleForLipidDoubleBondRowName(doubleBondInfo[0]);
//              this.drawARectangle(g2);
//            }
//          }
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
        //if click happened within heatmap boundaries
        if ((y>=imagePositionY_) && (y<imagePositionY_+renderedImage_.getHeight()) && (x>=imagePositionX_) && (x<(imagePositionX_+renderedImage_.getWidth()))) {
          int xInImage = x-imagePositionX_;
          int yInImage = y-imagePositionY_;
          if (heatmap_.getExpressionImageXStart()<=xInImage&&xInImage<heatmap_.getExpressionImageXEnd()&&
              heatmap_.getExpressionImageYStart()<=yInImage&&yInImage<heatmap_.getExpressionImageYEnd())
          {
          	int[] cellPos = heatmap_.getCellPosition(xInImage, yInImage);
            String experiment = heatmap_.getExperimentName(cellPos[0]);
            String analyte = heatmap_.getMolecularSpeciesLevelName(cellPos[1]);
            if (experiment!=null && analyte!=null)
            {
              ResultCompVO compVO = heatmap_.getCompVO(cellPos[0],cellPos[1]);
              if (compVO.getOriginalArea(0)>0 && heatMapListener_!=null){
                if (!isGrouped_){
                  if (e.getButton()==MouseEvent.BUTTON1){
                    if (!heatMapListener_.heatMapClicked(experiment,compVO,analyte,heatmap_.isMolecularSpeciesLevel(cellPos[1])));
                      this.mouseMoved(e);
                  }else if (e.getButton()==MouseEvent.BUTTON3 || e.isPopupTrigger()){
                  	Color attentionProbe = heatmap_.getAttentionProbe(cellPos[1], experiment);
                  	if (attentionProbe != null && attentionProbe == LipidomicsHeatMap.ATTENTION_COLOR_DOUBLE_PEAK || attentionProbe == LipidomicsHeatMap.ATTENTION_DOUBLE_AND_NOT_ALL_MODS)
                      this.mouseMoved(e);
                    else{
                      if (isAnalyteSelected(heatmap_.getAnalyteName(cellPos[1]),experimentNames_.get(cellPos[0]))){
                        selectSingleItem_.setEnabled(false);
                        deselectSingleItem_.setEnabled(true);
                      }else{
                        selectSingleItem_.setEnabled(true);
                        deselectSingleItem_.setEnabled(false);                  
                      }      
                      applySettingsPopup_.show(e.getComponent(), e.getX(), e.getY());
                      lastClickedCellPos_ = cellPos;
                    }  
                  }
                }else
                  this.mouseMoved(e);
              }else
                this.mouseMoved(e);
            }
          }else if (heatmap_.getRowNameStart()<xInImage&& xInImage<heatmap_.getRowNameEnd() &&
              heatmap_.getExpressionImageYStart()<=yInImage&&yInImage<heatmap_.getExpressionImageYEnd()){
            
          	int rowNumber = heatmap_.getRowNumber(x,y);
            String analyteName = heatmap_.getOriginalAnalyteName(rowNumber);
            boolean returnValue = false;
            if (isGrouped_)
              returnValue = heatMapListener_.analyteGroupClicked(analyteName,groupName_,maxIsotopes,rtTolerance_!=null,
              settingsVO_,heatmap_.getPreferredUnit(rowNumber),StaticUtils.getCorrespondingUnit(settingsVO_, heatmap_.getPreferredUnit(rowNumber),true));
            else{
              if (SwingUtilities.isLeftMouseButton(e)){
                returnValue = heatMapListener_.analyteClicked(analyteName,groupName_,maxIsotopes,rtTolerance_!=null,
                    settingsVO_,heatmap_.getPreferredUnit(rowNumber),StaticUtils.getCorrespondingUnit(settingsVO_, heatmap_.getPreferredUnit(rowNumber),true));
              } else if (SwingUtilities.isRightMouseButton(e)){
                if (isAnalyteSelected(rowNumber)){
                  selectItem_.setEnabled(false);
                  deselectItem_.setEnabled(true);
                }else{
                  selectItem_.setEnabled(true);
                  deselectItem_.setEnabled(false);                  
                }                  
                removeAnalytePopup_.show(e.getComponent(), e.getX(), e.getY());
                lastClickedRow_ = rowNumber;
              }
            } 
            if (!returnValue)
              this.mouseMoved(e);
          } else if (heatmap_.getExpressionImageXStart()<=xInImage&&xInImage<heatmap_.getExpressionImageXEnd()&&
              heatmap_.getColumnNameStart()<=yInImage&&yInImage<heatmap_.getColumnNameEnd()){
            String expName = heatmap_.getColumnName(x,y);
            if (e.getButton()==MouseEvent.BUTTON1){
              boolean returnValue = false;
              String preferredUnit = heatmap_.extractPreferredUnitForExp();
              if (isGrouped_)
                returnValue = heatMapListener_.experimentGroupClicked(expName,groupName_,maxIsotopes,rtTolerance_!=null,
                settingsVO_,preferredUnit,StaticUtils.getCorrespondingUnit(settingsVO_, preferredUnit,true));
              else
                returnValue = heatMapListener_.experimentClicked(expName,groupName_,maxIsotopes,rtTolerance_!=null,
                settingsVO_,preferredUnit,StaticUtils.getCorrespondingUnit(settingsVO_, preferredUnit,true));
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

  public void deselectAnalyte(int row)
  {
  	selectedMoleculeRows_.remove(row);
  }

  public boolean isAnalyteSelected(int row)
  {
    return selectedMoleculeRows_.contains(row);
  }

  public void selectAnalyte(int row)
  {
  	selectedMoleculeRows_.add(row);
  }

  /**
   * returns true when a single analyte was selected previously
   * @param analyteName the name of the analyte
   * @param expName the name of the experiment
   * @return true when a single analyte was selected previously
   */
  public boolean isAnalyteSelected(String analyteName, String expName)
  {
    if (selectedSingleMolecules_.containsKey(analyteName) && selectedSingleMolecules_.get(analyteName).containsKey(expName))
      return true;
    else
      return false;
  }

  /**
   * selects a single analyte in the heat map and stores additional information about this analyte
   * @param analyteName the name of the analyte
   * @param expName the name of the experiment
   * @param color the color of the expression image
   * @param modHash a hash table containing the possible modifications where an automated action can be executed
   * @param absolutePath the absolute file path to the experiment of the analte
   * @param updateableAndAnalyteBefore information about the experiments where the automated quantitation should take place
   * @param maxIsotope the highest possible isotope number for this analyte
   * @return true when the selection worked; false when the analyte was already selected
   */
  public boolean selectAnalyte(String analyteName, String expName, Color color, Hashtable<String,String> modHash, String absolutePath,
      Vector<AutoAnalyteAddVO> updateableAndAnalyteBefore, int maxIsotope)
  {
    if (selectedSingleMolecules_.containsKey(analyteName)){
      new WarningMessage(new JFrame(), "Error", "It is not allowed to select more than one from the same species");
      return false;
    }
    selectedSingleMolecules_.put(analyteName, new Hashtable<String,Color>());
    selectedSingleMoleculesMods_.put(analyteName, modHash);
    selectedSingleMolecules_.get(analyteName).put(expName, color);
    selectedSingleMoleculesAbsPaths_.put(analyteName, absolutePath);
    selectedSingleMoleculesMaxIso_.put(analyteName, maxIsotope);
    selectedSingleMoleculesAutoAnalyteAddVO_.put(analyteName, updateableAndAnalyteBefore);
    return true;
  }

  /**
   * selects a single analyte in the heat map and returns the previous color of this cell in the expression image
   * @param analyteName the name of the analyte
   * @param expName the name of the experiment
   * @return the previous color of this cell in the expression image
   */
  public Color deselectAnalyte(String analyteName, String expName)
  {
    Color color =  selectedSingleMolecules_.get(analyteName).get(expName);
    selectedSingleMolecules_.remove(analyteName);
    selectedSingleMoleculesMods_.remove(analyteName);
    selectedSingleMoleculesAbsPaths_.remove(analyteName);
    selectedSingleMoleculesAutoAnalyteAddVO_.remove(analyteName);
    selectedSingleMoleculesMaxIso_.remove(analyteName);  
    return color;
  }
  
  public int getSelectedIsotope(){
    int isotope = 0;
    if (maxIsotopes_.getSelectedItem()!=null)
      isotope = new Integer((String)maxIsotopes_.getSelectedItem()); 
    return isotope;
  }
  
  private ArrayList<String> getDisplayNames(){
  	ArrayList<String> names = new ArrayList<String>();
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
  
  public Hashtable<String,Hashtable<String,Vector<Double>>> getResultValues(String valueType) throws NumberFormatException, CalculationNotPossibleException{
    int exportType = ExportOptionsVO.EXPORT_NO_DEVIATION;
    if (isGrouped_) exportType = ExportOptionsVO.EXPORT_SD_DEV_AND_ERROR;
    ExportOptionsVO expOptions = new ExportOptionsVO(exportType,"1", false, true, false, false, 0, LipidomicsConstants.EXPORT_ANALYTE_TYPE_SPECIES);
    ResultDisplaySettingsVO settings = new ResultDisplaySettingsVO(settingsVO_);
    settings.setType(ResultDisplaySettingsVO.REL_VALUE);
    if (valueType!=null && valueType.length()>0)
      settings.setType(valueType);
    if (maxIsotopes_.getItemCount()>0)
      return HeatMapDrawing.extractValuesOfInterest(resultsOfOneGroup_, Integer.parseInt((String)maxIsotopes_.getSelectedItem()), settings, null, expOptions,modifications_);
    else
      return null;
  }
  
  
  /** returns the currently set export options, or the default settings if the
   * panel was not activated
   * @return the currently set export options, or the default settings if the panel was not activated
   */
  private ExportOptionsVO getExportOptions(){
    ExportOptionsVO expOptions = new ExportOptionsVO(ExportOptionsVO.EXPORT_NO_DEVIATION,null,true,false,false,false,6,LipidomicsConstants.EXPORT_ANALYTE_TYPE_SPECIES);
    if (exportSettings_ != null)
      expOptions = exportSettings_.getSettings();
    return expOptions;
  }
  
  
  /**
   * 
   * @return the type of export value used
   */
  public String getValueType(){
    return this.settingsVO_.getType();
  }
  
  protected boolean isMarkDoublePeaks()
  {
  	return markDoublePeaks_!=null && markDoublePeaks_.isSelected();
  }
  
  
  /**
   * returns a list of selected species plus additional information about this species
   * each piece of information is added as separate object to the return vector
   * first object:  sorted vector containing SimpleValueObject where the key is the analyte name (including RT), and the value is the experiment name
   * second object: a hash table containing all possible modifications for the selected analytes
   * third object:  sorted vector of absolute paths to the Excel result files of the selected single analyte (same order as first object)
   * fourth object: hash table, key: analyte name; value: vector of objects holding information about the other experiments where automated processing should take place
   * fifth object:  sorted vector of highest isotope numbers for each analyte (same order as first object)
   * @param currentlySelected SimpleValueObject of the currently selected species; label: analyte name; value: experiment name
   * @param modHash hash table containing the possible modifications where a quantitation can be performed; both the key and the value hold the same information: the modification
   * @param absolutePath the absolute paths to the Excel result files of the last clicked single analyte
   * @param updateableAndAnalyteBefore vector of objects holding information about the other experiments where automated processing should take place
   * @param maxIsotope highest isotope number for the currently selected analyte
   * @return a list of selected species plus additional information about this species where each piece of information is added as separate object to the return vector
   */
  @SuppressWarnings({ "rawtypes", "unchecked" })
  private Vector getSelectedIndividualSpeciesInCorrectOrder(SimpleValueObject currentlySelected, Hashtable<String,String> modHash,
      String absolutePath, Vector<AutoAnalyteAddVO> updateableAndAnalyteBefore, int maxIsotope){
    Vector result = new Vector();
    Vector<SimpleValueObject> sorted = new Vector<SimpleValueObject>();
    Hashtable<String,String> allMods = new Hashtable<String,String>();
    Vector<String> absPaths = new Vector<String>();
    Hashtable<String,Vector<AutoAnalyteAddVO>> updateablesAndAnalytesBefore = new Hashtable<String,Vector<AutoAnalyteAddVO>>();
    Vector<Integer> maxIsotopes = new Vector<Integer>();
    for (String molName :moleculeNames_){
      if (!selectedSingleMolecules_.containsKey(molName) && !currentlySelected.getLabel().equalsIgnoreCase(molName))
        continue;
      String expName;
      if (currentlySelected.getLabel().equalsIgnoreCase(molName)){
        expName = currentlySelected.getValue();
        for (String mod : modHash.keySet()) allMods.put(mod, mod);
        absPaths.add(absolutePath);
        updateablesAndAnalytesBefore.put(molName, updateableAndAnalyteBefore);
        maxIsotopes.add(maxIsotope);
      }else{
        expName = selectedSingleMolecules_.get(molName).keySet().iterator().next();
        for (String mod : selectedSingleMoleculesMods_.get(molName).keySet()) allMods.put(mod, mod);
        absPaths.add(selectedSingleMoleculesAbsPaths_.get(molName));
        updateablesAndAnalytesBefore.put(molName, selectedSingleMoleculesAutoAnalyteAddVO_.get(molName));
        maxIsotopes.add(selectedSingleMoleculesMaxIso_.get(molName));
      }
      sorted.add(new SimpleValueObject(molName,expName));
    }
    result.add(sorted);
    result.add(allMods);
    result.add(absPaths);
    result.add(updateablesAndAnalytesBefore);
    result.add(maxIsotopes);
    return result;
  }
  
  /**
   * Extracts LipidParameterSets belonging to the individual molecule names of one lipid class as displayed in the heat map
   * @param all all detected results of a certain experiment and lipid class
   * @param areaVO the value object corresponding to one data point (rectangle) in the heat map
   * @param msnOnly if true, only species verified by MSn evidence are exported
   * @return the original information relevant to this one data point (rectangle)
   */
  protected static Hashtable<String,Vector<LipidParameterSet>> getRelevantOriginalResults(Vector<LipidParameterSet> all, ResultAreaVO areaVO){
    Hashtable<String,Vector<LipidParameterSet>> results = new Hashtable<String,Vector<LipidParameterSet>>();
    for (LipidParameterSet set : all){
      if (!areaVO.getMoleculeNameWoRT().equalsIgnoreCase(set.getNameStringWithoutRt()))
        continue;
      if (!areaVO.belongsRtToThisAreaVO(set.getRt(), set.getModificationName()))
        continue;
      Vector<LipidParameterSet> ofOneMod = new Vector<LipidParameterSet>();
      if (results.containsKey(set.getModificationName()))
        ofOneMod = results.get(set.getModificationName());
      ofOneMod.add(set);
      results.put(set.getModificationName(), ofOneMod);
    }
    return results;
  }
  
  
  public void cleanup(){
  	ungroupedPartner_ = null;
    if (displaySettings_ != null)
    {
    	displaySettings_.removeActionListener(this);
    	displaySettings_.dispose();
      displaySettings_ = null;
    }
    if (selectionSettings_ != null) 
    {
    	selectionSettings_.removeActionListener(this);
    	selectionSettings_.dispose();
      selectionSettings_ = null;
    }
    if (combinedChartSettings_ != null)
    {
    	combinedChartSettings_.removeActionListener(this);
    	combinedChartSettings_.dispose();
      combinedChartSettings_ = null;
    }
    if (exportSettings_ != null) 
    {
    	exportSettings_.cleanup();
    	exportSettings_.dispose();
      exportSettings_ = null;
    }
    if (chromExport_ != null) 
    {
    	chromExport_.cleanup();
    	chromExport_.dispose();
      chromExport_ = null;
    }
    if (timer_!=null){
      timer_.cancel();
      timer_.purge();
    }
    timer_ = null;
  }
  
}
