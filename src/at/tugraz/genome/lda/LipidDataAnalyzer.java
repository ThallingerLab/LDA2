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

package at.tugraz.genome.lda;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.Toolkit;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.StringSelection;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Properties;
import java.util.Set;
import java.util.Timer;
import java.util.TimerTask;
import java.util.Vector;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import javax.imageio.ImageIO;
import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.ImageIcon;
import javax.swing.JApplet;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.ListSelectionModel;
import javax.swing.SwingConstants;
import javax.swing.UIManager;
import javax.swing.border.Border;
import javax.swing.border.EtchedBorder;
import javax.swing.border.TitledBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.TransformerConfigurationException;
import javax.xml.transform.TransformerException;

import uk.ac.ebi.pride.jmztab2.utils.errors.MZTabErrorType.Level;

import org.apache.batik.dom.GenericDOMImplementation;
import org.apache.batik.svggen.SVGGraphics2D;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;
import org.w3c.dom.DOMImplementation;
import org.w3c.dom.Document;

import at.tugraz.genome.lda.alex123.RdbOutputWriter;
import at.tugraz.genome.lda.analysis.AnalyteAddRemoveListener;
import at.tugraz.genome.lda.analysis.ClassNamesExtractor;
import at.tugraz.genome.lda.analysis.ComparativeAnalysis;
import at.tugraz.genome.lda.analysis.ComparativeNameExtractor;
import at.tugraz.genome.lda.analysis.ComparativeResultsLookup;
import at.tugraz.genome.lda.analysis.HeatMapClickListener;
import at.tugraz.genome.lda.analysis.exception.CalculationNotPossibleException;
import at.tugraz.genome.lda.exception.AbsoluteSettingsInputException;
import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.exception.ExcelInputFileException;
import at.tugraz.genome.lda.exception.ExportException;
import at.tugraz.genome.lda.exception.HydroxylationEncodingException;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.exception.NoRuleException;
import at.tugraz.genome.lda.exception.QuantificationException;
import at.tugraz.genome.lda.exception.RdbWriterException;
import at.tugraz.genome.lda.exception.RetentionTimeGroupingException;
import at.tugraz.genome.lda.exception.RulesException;
import at.tugraz.genome.lda.exception.SettingsException;
import at.tugraz.genome.lda.export.ExcelAndTextExporter;
import at.tugraz.genome.lda.export.LDAExporter;
import at.tugraz.genome.lda.export.OmegaCollector;
import at.tugraz.genome.lda.export.QuantificationResultExporter;
import at.tugraz.genome.lda.fragai.ExcelTargetListParser;
import at.tugraz.genome.lda.fragai.SpectraIdentifier;
import at.tugraz.genome.lda.fragai.SpectraInterpreter;
import at.tugraz.genome.lda.fragai.SpectraTextExporter;
import at.tugraz.genome.lda.fragai.SpectrumContainer;
import at.tugraz.genome.lda.glyco.FragGLiPanel;
import at.tugraz.genome.lda.fragai.CombinedSpectrumContainer;
import at.tugraz.genome.lda.interfaces.ColorChangeListener;
import at.tugraz.genome.lda.listeners.AnnotationThresholdListener;
import at.tugraz.genome.lda.masslist.MassListCreatorPanel;
import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
import at.tugraz.genome.lda.msn.MSnAnalyzer;
import at.tugraz.genome.lda.msn.RulesContainer;
import at.tugraz.genome.lda.msn.hydroxy.parser.HydroxyEncoding;
import at.tugraz.genome.lda.msn.vos.FattyAcidVO;
import at.tugraz.genome.lda.mztab.MztabUtils;
import at.tugraz.genome.lda.mztab.SmallMztabMolecule;
import at.tugraz.genome.lda.parser.LDAResultReader;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.quantification.LipidomicsAnalyzer;
import at.tugraz.genome.lda.quantification.QuantificationResult;
import at.tugraz.genome.lda.swing.AbsoluteQuantSettingsPanel;
import at.tugraz.genome.lda.swing.BarChartPainter;
import at.tugraz.genome.lda.swing.ClassesOverviewPanel;
import at.tugraz.genome.lda.swing.ColorChooserDialog;
import at.tugraz.genome.lda.swing.CutoffSettingsPanel;
import at.tugraz.genome.lda.swing.EditOmegaAssignmentJTable;
import at.tugraz.genome.lda.swing.EditRtDialog;
import at.tugraz.genome.lda.swing.ExportPanel;
import at.tugraz.genome.lda.swing.ExportSettingsPanel;
import at.tugraz.genome.lda.swing.GroupsPanel;
import at.tugraz.genome.lda.swing.HeatMapDrawing;
import at.tugraz.genome.lda.swing.InputDialog;
import at.tugraz.genome.lda.swing.JHyperlink;
import at.tugraz.genome.lda.swing.LipidomicsJTable;
import at.tugraz.genome.lda.swing.LipidomicsTableCellRenderer;
import at.tugraz.genome.lda.swing.LipidomicsTableModel;
import at.tugraz.genome.lda.swing.RangeColor;
import at.tugraz.genome.lda.swing.RecalculateMSnDialog;
import at.tugraz.genome.lda.swing.ResultDisplaySettings;
import at.tugraz.genome.lda.swing.ResultSelectionSettings;
import at.tugraz.genome.lda.swing.RuleDefinitionInterface;
import at.tugraz.genome.lda.swing.SpectrumUpdateListener;
import at.tugraz.genome.lda.target.JTargetFileWizard;
import at.tugraz.genome.lda.utils.ExcelUtils;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.verifier.DoubleVerifier;
import at.tugraz.genome.lda.vos.AbsoluteSettingsVO;
import at.tugraz.genome.lda.vos.AddAnalyteVO;
import at.tugraz.genome.lda.vos.AutoAnalyteAddVO;
import at.tugraz.genome.lda.vos.IntegerStringVO;
import at.tugraz.genome.lda.vos.QuantVO;
import at.tugraz.genome.lda.vos.RawQuantificationPairVO;
import at.tugraz.genome.lda.vos.ResultCompVO;
import at.tugraz.genome.lda.vos.ResultDisplaySettingsVO;
import at.tugraz.genome.lda.vos.ResultFileVO;
import at.tugraz.genome.lda.xml.AbsoluteQuantSettingsWholeReader;
import at.tugraz.genome.lda.xml.AbsoluteQuantSettingsWholeWriter;
import at.tugraz.genome.lda.xml.AbstractXMLSpectraReader;
import at.tugraz.genome.lda.xml.CutoffSettingsReader;
import at.tugraz.genome.lda.xml.CutoffSettingsWriter;
import at.tugraz.genome.lda.xml.RawToChromThread;
import at.tugraz.genome.maspectras.GlobalConstants;
import at.tugraz.genome.maspectras.chromaviewer.MSMapViewer;
import at.tugraz.genome.maspectras.chromaviewer.MSMapViewerFactory;
import at.tugraz.genome.maspectras.parser.exceptions.SpectrummillParserException;
import at.tugraz.genome.maspectras.quantification.CgAreaStatus;
import at.tugraz.genome.lda.quantification.LipidomicsDefines;
import at.tugraz.genome.maspectras.quantification.CgException;
import at.tugraz.genome.maspectras.quantification.CgProbe;
import at.tugraz.genome.maspectras.quantification.ChromatogramHeaderFileReader;
import at.tugraz.genome.maspectras.quantification.ChromatogramReader;
import at.tugraz.genome.maspectras.utils.Calculator;
import at.tugraz.genome.maspectras.utils.StringUtils;
import at.tugraz.genome.util.index.IndexFileException;
import at.tugraz.genome.voutils.GeneralComparator;

import org.jogamp.java3d.utils.applet.MainFrame;

import de.isas.mztab2.io.MzTabValidatingWriter;
import de.isas.mztab2.io.MzTabWriter;
import de.isas.mztab2.io.MzTabWriterDefaults;
import de.isas.mztab2.model.Assay;
import de.isas.mztab2.model.CV;
import de.isas.mztab2.model.Contact;
import de.isas.mztab2.model.Database;
import de.isas.mztab2.model.Instrument;
import de.isas.mztab2.model.Metadata;
import de.isas.mztab2.model.MsRun;
import de.isas.mztab2.model.MzTab;
import de.isas.mztab2.model.Parameter;
import de.isas.mztab2.model.Sample;
import de.isas.mztab2.model.SampleProcessing;
import de.isas.mztab2.model.SmallMoleculeEvidence;
import de.isas.mztab2.model.SmallMoleculeFeature;
import de.isas.mztab2.model.SmallMoleculeSummary;
import de.isas.mztab2.model.Software;
import de.isas.mztab2.model.StudyVariable;
import de.isas.mztab2.model.ValidationMessage;

/**
 * 
 * @author Juergen Hartler
 * @author Leonida M. Lamp
 *
 */
public class LipidDataAnalyzer extends JApplet implements ActionListener,HeatMapClickListener,AnalyteAddRemoveListener,ColorChangeListener,SpectrumUpdateListener
{
    
  private static final long serialVersionUID = -1358973418774645198L;
  
  private static final int FRAME_WIDTH = 1050;
  private static final int FRAME_HEIGHT = 950;  
  private static final int RESULTS_TABLE_WIDTH = 950;
  
  private static MainFrame frame_;
  private JFileChooser mzXMLFileChooser;
  private JFileChooser mzXMLDirChooser;
  private JTabbedPane mainTabs;
  private QuantificationMenu singleQuantMenu_;
  private QuantificationMenu batchQuantMenu_;
  private JTabbedPane quantitationPane_;
  private JPanel resultsMenu_;
  private JPanel resultsPanel_;
  private JPanel settingsPanel_;
  private JPanel licensePanel_;
  private MassListCreatorPanel massListPanel_;
  private JTargetFileWizard targetFilePanel_;
  private FragGLiPanel fragGLiPanel_;
  
  /**
   * what do I want and need? I need to read in 
   */
  
  private JPanel helpPanel_;
  private JPanel aboutPanel_;
  private JPanel displayTopMenu;
  private QuantificationThread quantThread_;
  private BatchQuantThread batchQuantThread_;
  private RawToMzxmlThread rawmzThread_;
  private MzxmlToChromThread mzToChromThread_;
  private Timer timer_;
  
  private JLabel resultStatus_;

  private JTabbedPane resultTabs_;
  private JPanel resultStatusPanel_;
  private JPanel resultsSelectionPanel_;
  private Hashtable<String,JTabbedPane> molBarCharts_; 

  private JButton jButtonResultsDirOpen_;
  private JButton jButtonResultFilesOpen_;
  private JButton jButtonResultFilesRemove_;
  private JButton jButtonResultAddToGroup_;
  private JCheckBox separateHitsByRT_;
  private JTextField rtGroupingTime_;
  private JLabel rtTimeUnit_;
  private JButton jButtonResultFilesClean_;
  private GroupsPanel groupsPanel_;
  private AbsoluteQuantSettingsPanel quantSettingsPanel_;
  private CutoffSettingsPanel cutoffSettingsPanel_;
  private JButton jButtonResultFilesAccept_;
  private JButton jButtonResultAbsQuant_;
  private JButton jButtonResultCutoff_;
  private JFileChooser resultsDirChooser_;
  private JFileChooser resultFilesChooser_;
  private JScrollPane analysisTablePane_;
  private JPanel analysisSelectionTablePanel_;
  private JTable resultFilesDisplayTable_;
  private ListSelectionModel resultListSelectionModel_;
  private ImageIcon addFilesIcon_ = new ImageIcon(getClass().getResource(
  "/images/addFiles.gif"));
  private ImageIcon addFolderIcon_ = new ImageIcon(getClass().getResource(
  "/images/addFromFolder.gif"));
  private ImageIcon removeFilesIcon_ = new ImageIcon(getClass().getResource(
  "/images/removeFiles.gif"));
  private ImageIcon addToGroupIcon_ = new ImageIcon(getClass().getResource(
  "/images/User16.gif"));
  private ImageIcon ldaLogo_ = new ImageIcon(getClass().getResource(
      "/images/lda_logo.png"));
  private ImageIcon tugLogo_ = new ImageIcon(getClass().getResource(
      "/images/logo-RGB.png"));

  private JTextField internalStandardSelection_;
  private JTextField externalStandardSelection_;
  private JTextField correctOrderFile_;
  private JFileChooser correctOrderFileChooser_;
  
  private JTextField selectedChromFile;
  private JButton jButtonChromOpen;
  private JFileChooser chromFileChooser_;
  private JFileChooser quantFileChooser_;
  private JFileChooser quantDirChooser_;
  private JComboBox<String> ionModeOrder_;

  
  private JTextField selectedResultFile;
  private JButton jButtonResultOpen;
  private JButton startDisplay;
  private JFileChooser resultFileChooser_;
  private LipidomicsJTable displayTable_;
//  private BatchQuantificationTable batchQuantTable_;
//  private BatchQuantificationTableModel batchQuantTableModel_;
  private ListSelectionModel listSelectionModel;
  private JPanel selectionPane;
  private JPanel tableContainer;
  private JPanel tablePanel_;
  private DisplayTolerancePanel displayTolerancePanel_;
  private JComboBox<String> selectedSheet_;

  private Lipidomics2DPainter l2DPainter_;
  private Lipidomics2DPainter spectrumPainter_;
  private LipidomicsItemListener changeIsotopeListener_;
  
  private JPanel displayPanel_;
  private JSplitPane majorSplitPane_;
  private JSplitPane topSplitPane_;
  private RuleDefinitionInterface userInterface;
  private JScrollPane tablePane;
  private QuantificationResult result_;
  private QuantificationResult originalResult_;
  private QuantificationResult mSnResult_;
  private QuantificationResult chainResult_;
  private Hashtable<Integer,Integer> resultPositionToOriginalLoopkup_ = new Hashtable<Integer,Integer>();
  private Hashtable<Integer,String> resultPositionToMolecularSpeciesLookup_;
  private Hashtable<String,Integer> orderResultsType_ = new Hashtable<String,Integer>();
  private Hashtable<String,Boolean> resultsShowModification_ = new Hashtable<String,Boolean>();
  private Vector<File> resultFiles_;
  private int currentSelected_ = -1;
  private String currentSelectedSheet_ = "";
  private ChromatogramReader reader_;
  private LipidomicsAnalyzer analyzer_;
  private JLabel resultWarningLabel_;
  private JLabel cutoffWarningLabel_;
  private JButton resultLoadButton_; 
  private JButton resultSaveButton_; 
  private JButton cutoffLoadButton_; 
  private JButton cutoffSaveButton_; 

  private JFileChooser saveAbsSettingsFileChooser_;
  private JFileChooser loadAbsSettingsFileChooser_;

  private JFileChooser saveCutoffSettingsFileChooser_;
  private JFileChooser loadCutoffSettingsFileChooser_;
  
  private LipidParameterSet params_;
  private MSMapViewer viewer_;
  private JPanel l2dPanel_;
  /** layout for the 2D-panel - used to remove the Component at CENTER position*/
  private BorderLayout l2dPanelLayout_;
  private JPanel spectrumPanel_;
  private boolean readFromRaw_;
  
  private JRadioButton m_chkRaw_;
  private JComboBox<String> isotope_;
  private JLabel isotopeLabel_;
  private JRadioButton m_chkSmooth_;
  private JButton m_upButton_;
  private JButton m_dnButton_;
  private JLabel lx_min_;
  private JLabel lx_max_;
  private JTextField m_minTimeText_;
  private JTextField m_maxTimeText_;
  private JLabel lz_min_;
  private JLabel lz_max_;
  private JTextField mz_minTimeText_;
  private JTextField mz_maxTimeText_;
  private JButton m_zoomIn_;
  private JButton m_zoomAll_;
  private JButton mz_zoomIn_;
  private JButton mz_zoomAll_;
  private JRadioButton relAbund_;
  private JRadioButton absAbund_;
  private JLabel spectrumSelectedLabel_;
  private JLabel spectrumSelected_;
  private JLabel rtSelectedLabel_;
  private JLabel rtSelected_;
  private JLabel precursorSelectedLabel_;
  private JLabel[] precursorSelected_;
  private JLabel msLevelSelectedLabel_;
  private JLabel msLevelSelected_;
  private JButton spectrumEarlier_;
  private JButton spectrumLater_;
  private JLabel annotationLabel_;
  private JTextField annotationThreshold_;
  private JLabel annotationUnit_;
  private ExportPanel exportSpectra_;
  private JFileChooser exportFileChooser_;
  
  private JButton storeSelectedAreas_;
  
  private ClassesOverviewPanel classOverviewPanel_;
  private ClassesOverviewPanel classOverviewGroupPanel_;
  
  private ComparativeAnalysis analysisModule_;
  
  private Hashtable<String,HeatMapDrawing> heatmaps_;
  private Hashtable<String,HeatMapDrawing> groupHeatmaps_;
  private Hashtable<String,String> expDisplayNamesLookup_;
  private Hashtable<String,String> groupDisplayNamesLookup_;
  
  private ColorChooserDialog colorChooserDialog_;
  
  private boolean displaysMs2_;
//  private int currentMs2Position_;
  private Vector<Vector<CgProbe>> ms1ProbesWhileMs2Display_;
  
  /** panel for the selection of fragmentation language settings */
  private JPanel inputFragSettings2_;
  /** combo box containing the available MS machines*/
  private JComboBox<String> msMachineTypes_;
  /** combo box containing the fragmentation settings for each machine */
  private JComboBox<String> fragmentationSettings1_;
  /** combo box containing the fragmentation settings for each machine for a second selection*/
  private JComboBox<String> fragmentationSettings2_;
  
  private final static String CHANGE_SEPARATE_RT_STATUS = "rtActivateTextBox";
  
  /** Global finalButtonSection to delete it in the makeDisplayRemoveOperations */
  private JPanel finalButtonSection_;
  
  /** Object for the Rule Definition Interface */
  private RuleDefinitionInterface msnUserInterfaceObject_;
  
  private final static String DEFAULT_ANNOTATION_CUTOFF = "5";
  
  protected final static Font SELECT_FIELD_FONT = new Font("Helvetica",Font.PLAIN,10);
  
  /** a dialog field showing the new MSn assignment*/
  private RecalculateMSnDialog recalcDialog_;
  
  /** a dialog field for changing the retention time of a hit*/
  private EditRtDialog editRtDialog_;
  /** true for the display of shotgun data*/
  private boolean shotgunIsDisplayed_;
  
  /** joined panel for exporting things from heat maps*/
  private ExportSettingsPanel exportSettings_ = null;
  /** joined panel for exporting things from group heat maps*/
  private ExportSettingsPanel exportSettingsGroup_ = null;
  /** the lock range must be updated if true*/
  private boolean lockRangeUpdateRequired_ = false;
  
  /**for oxLipids */
  private JCheckBox showMSnEvidenceStat_;
  private JCheckBox showChainEvidenceStat_;
  private JCheckBox combineWithOx_;
  private boolean combineOxWithNonOx_;
  
  private int statisticsViewMode_ = 0;
  //TODO: set to false for production version
  private boolean exportChromatogramsFromDRView_ = false;
  
  public LipidDataAnalyzer(){
  	this.checkLicense();
    this.createDisplayTopMenu();
    this.batchQuantMenu_ = new QuantificationMenu(true,this);
    this.singleQuantMenu_ = new QuantificationMenu(false,this);
    this.singleQuantMenu_.setBatchQuantTableForSingleQuant(batchQuantMenu_);
    this.createResultsMenu();
    this.initL2dPanel();
    
    displaysMs2_ = false;
    shotgunIsDisplayed_ = false;
    
    selectionPane = new JPanel();
    selectionPane.setLayout(new BoxLayout(selectionPane, BoxLayout.LINE_AXIS));
    Vector<LipidParameterSet> dummy = new Vector<LipidParameterSet>();
    
    displayTable_ = new LipidomicsJTable(new LipidomicsTableModel(dummy,dummy,false,false,false),new LipidomicsTableCellRenderer(),false,0,true,this);
    listSelectionModel = displayTable_.getSelectionModel();
    displayTable_.setSelectionModel(listSelectionModel);
    listSelectionModel.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
    
    JPanel listContainer = new JPanel(new GridLayout(1,1));
    tablePane = new JScrollPane(displayTable_);

    displayTolerancePanel_ = new DisplayTolerancePanel(this,true);
    
    tablePanel_ = new JPanel();
    tablePanel_.setLayout(new BorderLayout());
    selectedSheet_ = new JComboBox<String>();
    tablePanel_.add(selectedSheet_,BorderLayout.NORTH);
    tablePanel_.add(tablePane,BorderLayout.CENTER);

    
    tableContainer = new JPanel(new GridLayout(1,1));
    tableContainer.setBorder(BorderFactory.createTitledBorder("Results"));
    tableContainer.add(tablePanel_);
    tablePane.setPreferredSize(new Dimension(500, 130));  
    
    selectionPane.add(tableContainer);
    selectionPane.add(listContainer);
    

    selectionPane.setMinimumSize(new Dimension(270, 50));
    selectionPane.setPreferredSize(new Dimension(240, 110));
    selectionPane.setVisible(false);
    
    mainTabs = new JTabbedPane();
    
    
    
    displayPanel_ = new JPanel();
    displayPanel_.setLayout(new BorderLayout());
    displayPanel_.add(displayTopMenu,BorderLayout.NORTH);
    displayPanel_.add(selectionPane,BorderLayout.WEST);

    JPanel singleQuantificationPanel = new JPanel();
    singleQuantificationPanel.setLayout(new BorderLayout());
    JPanel batchQuantificationPanel = new JPanel();
    batchQuantificationPanel.setLayout(new BorderLayout());
    singleQuantificationPanel.add(singleQuantMenu_);
    batchQuantificationPanel.add(batchQuantMenu_);
    resultsPanel_ = new JPanel();
    resultsPanel_.setLayout(new BorderLayout());
    settingsPanel_ = new JPanel();
    initSettingsPanel();
    licensePanel_ = new JPanel();
    massListPanel_ = new MassListCreatorPanel();
    targetFilePanel_ = new JTargetFileWizard();
    fragGLiPanel_ = new FragGLiPanel();
    helpPanel_ = new JPanel();
    initHelpPanel();
    aboutPanel_ = new JPanel();
    initAboutPanel();
    
    quantitationPane_ = new JTabbedPane();
    quantitationPane_.setOpaque(true);
    quantitationPane_.setBackground(Color.decode("#EEEEEE"));
    quantitationPane_.addTab("Batch Quantitation", batchQuantificationPanel);
    quantitationPane_.setToolTipTextAt(quantitationPane_.indexOfComponent(batchQuantificationPanel), TooltipTexts.TABS_MAIN_BATCH);
    quantitationPane_.addTab("Single Quantitation", singleQuantificationPanel);
    quantitationPane_.setToolTipTextAt(quantitationPane_.indexOfComponent(singleQuantificationPanel), TooltipTexts.TABS_MAIN_SINGLE);
    mainTabs.addTab("Quantitation", quantitationPane_);
    mainTabs.setToolTipTextAt(mainTabs.indexOfComponent(quantitationPane_), TooltipTexts.TABS_MAIN_QUANTITATION);
    mainTabs.addTab("Statistical Analysis", resultsPanel_);
    mainTabs.setToolTipTextAt(mainTabs.indexOfComponent(resultsPanel_), TooltipTexts.TABS_MAIN_STATISTICS);
    mainTabs.addTab("Display Results", displayPanel_);
    mainTabs.setToolTipTextAt(mainTabs.indexOfComponent(displayPanel_), TooltipTexts.TABS_MAIN_DISPLAY);
    if (Settings.SHOW_LCCL)
    {
    	mainTabs.addTab("MassList Creator", massListPanel_);
      mainTabs.setToolTipTextAt(mainTabs.indexOfComponent(massListPanel_), TooltipTexts.TABS_MAIN_MASSLIST_CREATOR);
    	mainTabs.addTab("LC=CL", targetFilePanel_);
      mainTabs.setToolTipTextAt(mainTabs.indexOfComponent(targetFilePanel_), TooltipTexts.TABS_MAIN_TARGET);
    }
    if (Settings.SHOW_FRAGGLI)
    {
    	mainTabs.addTab("FraGLi", fragGLiPanel_);
    }
    mainTabs.addTab("Settings", settingsPanel_);
    mainTabs.setToolTipTextAt(mainTabs.indexOfComponent(settingsPanel_), TooltipTexts.TABS_MAIN_SETTINGS);
    if (LicenseChecker.isCheckLicense())
    {
      mainTabs.addTab("License", licensePanel_);
      mainTabs.setToolTipTextAt(mainTabs.indexOfComponent(licensePanel_), TooltipTexts.TABS_MAIN_LICENSE);
      LicenseChangeListener licenseListener = new LicenseChangeListener();
      mainTabs.addChangeListener(licenseListener);
    }
    mainTabs.addTab("Help", helpPanel_);
    mainTabs.setToolTipTextAt(mainTabs.indexOfComponent(helpPanel_), TooltipTexts.TABS_MAIN_HELP);
    mainTabs.addTab("About", aboutPanel_);
    mainTabs.setToolTipTextAt(mainTabs.indexOfComponent(aboutPanel_), TooltipTexts.TABS_MAIN_ABOUT);
    mainTabs.setSelectedIndex(0);
    resultTabs_= new JTabbedPane();
    resultsPanel_.add(resultTabs_,BorderLayout.CENTER);
    resultStatusPanel_ = new JPanel(new BorderLayout());
    resultStatus_ = new JLabel();
    resultStatusPanel_.add(resultStatus_,BorderLayout.WEST);
    resultsPanel_.add(resultStatusPanel_,BorderLayout.SOUTH);
    
    resultsSelectionPanel_ = new JPanel();
    resultsSelectionPanel_.setLayout(new BorderLayout());
    JScrollPane scrollPane = new JScrollPane(resultsMenu_);
    resultsSelectionPanel_.add(scrollPane);
    resultTabs_.addTab("Selection",resultsSelectionPanel_);
    resultTabs_.setToolTipTextAt(0, TooltipTexts.TABS_RESULTS_SELECTION);
//    resultsPanel_.add(resultsMenu_);
    this.add(mainTabs);
    this.initTimer();
  }
  
  private void initL2dPanel(){
    l2dPanel_ = new JPanel();
    l2dPanelLayout_ = new BorderLayout();
    l2dPanel_.setLayout(l2dPanelLayout_);
    JPanel l2dMenu = new JPanel();
    l2dMenu.setPreferredSize(new Dimension(95,50));
    l2dPanel_.add(l2dMenu,BorderLayout.EAST);
    l2dMenu.setLayout(new GridBagLayout());
    Font smallFont=new Font("Helvetica",Font.PLAIN,9);
    Font toolBarFont=new Font("Helvetica",Font.PLAIN,10);
    
//    JPanel isotopePanel = new JPanel();
//    isotopePanel.setLayout(new BorderLayout());
    isotopeLabel_ = new JLabel("Isotope:");
    isotopeLabel_.setFont(smallFont);
    l2dMenu.add(isotopeLabel_,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
      ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(1, 6, 1, 0), 0, 0));
//    isotopePanel.add(isotopeLabel,BorderLayout.EAST);
    isotope_ = new JComboBox<String>();
//    isotope_.addItem("0");
//    isotope_.addItem("1");
//    isotope_.addItem("2");
//    isotope_.addItem("3");
    isotope_.setFont(smallFont);
    changeIsotopeListener_ = new LipidomicsItemListener("ChangeIsotope");
    isotope_.addItemListener(changeIsotopeListener_);
    isotope_.setSize(isotope_.getWidth()*4, isotope_.getHeight());
    l2dMenu.add(isotope_,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(1, 0, 1, 0), 0, 0));
//     isotopePanel.add(isotope_,BorderLayout.CENTER);
//    l2dMenu.add(isotopePanel,new GridBagConstraints(0, 0, 1, 2, 0.0, 0.0
//        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(1, 6, 1, 0), 0, 0));
    
    
    ButtonGroup group = new ButtonGroup();
    m_chkRaw_ = new JRadioButton("Raw");
    group.add(m_chkRaw_);
    m_chkRaw_.setFont(smallFont);
    m_chkRaw_.addItemListener(new LipidomicsItemListener("DisplayModeRaw"));
    l2dMenu.add(m_chkRaw_,new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(1, 2, 1, 0), 0, 0));
    m_chkSmooth_ = new JRadioButton("Smooth");
    m_chkSmooth_.setFont(smallFont);
    m_chkSmooth_.setSelected(true);
    m_chkSmooth_.addItemListener(new LipidomicsItemListener("DisplayModeSmooth"));
    group.add(m_chkSmooth_);
    l2dMenu.add(m_chkSmooth_,new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(1, 2, 1, 0), 0, 0));

    m_upButton_ = new JButton("- mz");
    m_upButton_.setMargin(new Insets(1,1,1,2));
    m_upButton_.setFont(smallFont);
    m_upButton_.setActionCommand("Dn");
    m_upButton_.addActionListener(this);
    l2dMenu.add(m_upButton_,new GridBagConstraints(1, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(1, 0, 1, 0), 2, 1));
      
    
    m_dnButton_ = new JButton("+ mz");
    m_dnButton_.setMargin(new Insets(1,1,1,1));
    m_dnButton_.setFont(smallFont);
    m_dnButton_.setActionCommand("Up");
    m_dnButton_.addActionListener(this);
    l2dMenu.add(m_dnButton_,new GridBagConstraints(1, 2, 1, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(1, 0, 1, 0), 2, 1));    

    // **** Zoom ****
    JPanel l2dZoomPanel = new JPanel();
    l2dZoomPanel.setLayout(new GridBagLayout());

    l2dMenu.add(l2dZoomPanel,new GridBagConstraints(0, 7, 2, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(1, 0, 0, 0), 0, 0));

    lx_min_ = new JLabel();
    lx_min_.setText("t[min]:");
    lx_min_.setFont(toolBarFont);
    l2dZoomPanel.add(lx_min_,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(1, 3, 1, 0), 0, 0));
    m_minTimeText_ = new JTextField(3);
    m_minTimeText_.setText("");
    m_minTimeText_.setFont(toolBarFont);
    m_minTimeText_.setHorizontalAlignment(JTextField.RIGHT);
    l2dZoomPanel.add(m_minTimeText_,new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    lx_max_ = new JLabel();
    lx_max_.setText("t[max]:");
    lx_max_.setFont(toolBarFont);
    l2dZoomPanel.add(lx_max_,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(1, 3, 1, 0), 0, 0));
    m_maxTimeText_ = new JTextField(3);
    m_maxTimeText_.setText("");
    m_maxTimeText_.setFont(toolBarFont);
    m_maxTimeText_.setHorizontalAlignment(JTextField.RIGHT);
    l2dZoomPanel.add(m_maxTimeText_,new GridBagConstraints(1, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    
    m_zoomIn_ = new JButton("Zoom in");
    m_zoomIn_.setMargin(new Insets(1,5,1,5));
    m_zoomIn_.setFont(toolBarFont);
    m_zoomIn_.setActionCommand("ZoomIn");
    m_zoomIn_.addActionListener(this);
    l2dMenu.add(m_zoomIn_,new GridBagConstraints(0, 8, 2, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(3, 0, 0, 0), 2, 0));
//    m_comPnl.add(m_zoomIn,new GridBagConstraints(0, 14, 2, 1, 0.0, 0.0
//        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    m_zoomAll_ = new JButton("Zoom all");
    m_zoomAll_.setMargin(new Insets(1,5,1,5));
    m_zoomAll_.setFont(toolBarFont);
    m_zoomAll_.setActionCommand("ZoomAll");
    m_zoomAll_.addActionListener(this);
    l2dMenu.add(m_zoomAll_,new GridBagConstraints(0, 9, 2, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(3, 0, 0, 0), 0, 0));

    storeSelectedAreas_ = new JButton("Save");
    storeSelectedAreas_.setMargin(new Insets(1,5,1,5));
    storeSelectedAreas_.setFont(toolBarFont);
    storeSelectedAreas_.setActionCommand("SaveSelectedAreas");
    storeSelectedAreas_.addActionListener(this);
    l2dMenu.add(storeSelectedAreas_,new GridBagConstraints(0, 10, 2, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(15, 0, 0, 0), 0, 0));

    spectrumPanel_ = new JPanel();
    spectrumPanel_.setLayout(new BorderLayout());
    JPanel spectrumMenu = new JPanel();
    spectrumMenu.setPreferredSize(new Dimension(95,80));
    spectrumPanel_.add(spectrumMenu,BorderLayout.EAST);
    spectrumMenu.setLayout(new GridBagLayout());
    
    JPanel l2dSpectrumSelectionPanel = new JPanel();
    spectrumMenu.add(l2dSpectrumSelectionPanel,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(1, 0, 0, 0), 0, 0));
    l2dSpectrumSelectionPanel.setLayout(new GridBagLayout());

    spectrumSelectedLabel_ = new JLabel();
    spectrumSelectedLabel_.setText("Spect.: ");
    spectrumSelectedLabel_.setFont(toolBarFont);
    l2dSpectrumSelectionPanel.add(spectrumSelectedLabel_,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(1, 3, 1, 0), 0, 0));
    spectrumSelected_ = new JLabel();
    spectrumSelected_.setText("");
    spectrumSelected_.setFont(toolBarFont);
    l2dSpectrumSelectionPanel.add(spectrumSelected_,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(1, 3, 1, 0), 0, 0));
    rtSelectedLabel_ = new JLabel();
    rtSelectedLabel_.setText("RT: ");
    rtSelectedLabel_.setFont(toolBarFont);
    l2dSpectrumSelectionPanel.add(rtSelectedLabel_,new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(1, 3, 1, 0), 0, 0));
    rtSelected_ = new JLabel();
    rtSelected_.setText("");
    rtSelected_.setFont(toolBarFont);
    l2dSpectrumSelectionPanel.add(rtSelected_,new GridBagConstraints(1, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(1, 3, 1, 0), 0, 0));

    precursorSelectedLabel_ = new JLabel();
    precursorSelectedLabel_.setText("Prec: ");
    precursorSelectedLabel_.setFont(toolBarFont);
    l2dSpectrumSelectionPanel.add(precursorSelectedLabel_,new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(1, 3, 1, 0), 0, 0));
    precursorSelected_ = new JLabel[5];
    for (int i=0; i!=precursorSelected_.length;i++){
      precursorSelected_[i] = new JLabel();
      precursorSelected_[i].setText("");
      precursorSelected_[i].setFont(toolBarFont);
      l2dSpectrumSelectionPanel.add(precursorSelected_[i],new GridBagConstraints(1, 2+i, 1, 1, 0.0, 0.0
          ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(1, 3, 1, 0), 0, 0));
    }

    msLevelSelectedLabel_ = new JLabel();
    msLevelSelectedLabel_.setText("Level: ");
    msLevelSelectedLabel_.setFont(toolBarFont);
    l2dSpectrumSelectionPanel.add(msLevelSelectedLabel_,new GridBagConstraints(0, 7, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(1, 3, 1, 0), 0, 0));
    msLevelSelected_ = new JLabel();
    msLevelSelected_.setText("");
    msLevelSelected_.setFont(toolBarFont);
    l2dSpectrumSelectionPanel.add(msLevelSelected_,new GridBagConstraints(1, 7, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(1, 3, 1, 0), 0, 0));
    spectrumEarlier_ = new JButton(" - ");
    spectrumEarlier_.setMargin(new Insets(1,1,1,1));
    spectrumEarlier_.setFont(smallFont);
    spectrumEarlier_.setActionCommand("SpectMinus");
    spectrumEarlier_.addActionListener(this);
    l2dSpectrumSelectionPanel.add(spectrumEarlier_,new GridBagConstraints(0, 8, 1, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(1, 4, 1, 0), 0, 0));
    spectrumLater_ = new JButton(" + ");
    spectrumLater_.setMargin(new Insets(1,1,1,1));
    spectrumLater_.setFont(smallFont);
    spectrumLater_.setActionCommand("SpectPlus");
    spectrumLater_.addActionListener(this);
    l2dSpectrumSelectionPanel.add(spectrumLater_,new GridBagConstraints(1, 8, 1, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(1, 4, 1, 0), 0, 0));
    
    JPanel annotationPanel = new JPanel();
    annotationLabel_ = new JLabel("Annot.:");
    annotationLabel_.setFont(smallFont);
    annotationPanel.add(annotationLabel_);
    annotationThreshold_ = new JTextField(2);
    annotationThreshold_.setInputVerifier(new DoubleVerifier());
    annotationThreshold_.setFont(smallFont);
    annotationThreshold_.setText(DEFAULT_ANNOTATION_CUTOFF);
    annotationThreshold_.setHorizontalAlignment(JTextField.RIGHT);
    AnnotationThresholdListener annotationListener = new AnnotationThresholdListener(this);
    annotationThreshold_.getDocument().addDocumentListener(annotationListener); 
    annotationThreshold_.addFocusListener(annotationListener);     

    annotationPanel.add(annotationThreshold_);
    annotationUnit_ = new JLabel("%");
    annotationUnit_.setFont(smallFont);
    annotationPanel.add(annotationUnit_);
    
    l2dSpectrumSelectionPanel.add(annotationPanel,new GridBagConstraints(0, 9, 2, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(1, 3, 1, 0), 0, 0));
    
    
    JPanel l2dSpectrumIntensityPanel = new JPanel();
    spectrumMenu.add(l2dSpectrumIntensityPanel,new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(1, 0, 0, 0), 0, 0));
    l2dSpectrumIntensityPanel.setLayout(new GridBagLayout());
    group = new ButtonGroup();
    relAbund_ = new JRadioButton("rel. Abund.");
    group.add(relAbund_);
    relAbund_.setFont(smallFont);
    relAbund_.setSelected(true);
    relAbund_.addItemListener(new LipidomicsItemListener("DisplayModeAbundance"));
    l2dSpectrumIntensityPanel.add(relAbund_,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(1, 2, 1, 0), 0, 0));
    absAbund_ = new JRadioButton("abs. Abund.");
    absAbund_.setFont(smallFont);
    absAbund_.addItemListener(new LipidomicsItemListener("DisplayModeAbundance"));
    group.add(absAbund_);
    l2dSpectrumIntensityPanel.add(absAbund_,new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(1, 2, 1, 0), 0, 0));
    
    // **** Zoom ****
    JPanel spectrumZoomPanel = new JPanel();
    spectrumZoomPanel.setLayout(new GridBagLayout());
    spectrumMenu.add(spectrumZoomPanel,new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(1, 0, 0, 0), 0, 0));
    lz_min_ = new JLabel();
    lz_min_.setText("mz[min]:");
    lz_min_.setFont(toolBarFont);
    spectrumZoomPanel.add(lz_min_,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(1, 3, 1, 0), 0, 0));
    mz_minTimeText_ = new JTextField(3);
    mz_minTimeText_.setText("");
    mz_minTimeText_.setFont(toolBarFont);
    mz_minTimeText_.setHorizontalAlignment(JTextField.RIGHT);
    spectrumZoomPanel.add(mz_minTimeText_,new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    lz_max_ = new JLabel();
    lz_max_.setText("mz[max]:");
    lz_max_.setFont(toolBarFont);
    spectrumZoomPanel.add(lz_max_,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(1, 3, 1, 0), 0, 0));
    mz_maxTimeText_ = new JTextField(3);
    mz_maxTimeText_.setText("");
    mz_maxTimeText_.setFont(toolBarFont);
    mz_maxTimeText_.setHorizontalAlignment(JTextField.RIGHT);
    spectrumZoomPanel.add(mz_maxTimeText_,new GridBagConstraints(1, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    mz_zoomIn_ = new JButton("Zoom in");
    mz_zoomIn_.setMargin(new Insets(1,5,1,5));
    mz_zoomIn_.setFont(toolBarFont);
    mz_zoomIn_.setActionCommand("ZoomMzIn");
    mz_zoomIn_.addActionListener(this);
    spectrumZoomPanel.add(mz_zoomIn_,new GridBagConstraints(0, 2, 2, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(3, 0, 0, 0), 2, 0));
//    m_comPnl.add(m_zoomIn,new GridBagConstraints(0, 14, 2, 1, 0.0, 0.0
//        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    mz_zoomAll_ = new JButton("Zoom all");
    mz_zoomAll_.setMargin(new Insets(1,5,1,5));
    mz_zoomAll_.setFont(toolBarFont);
    mz_zoomAll_.setActionCommand("ZoomMzAll");
    mz_zoomAll_.addActionListener(this);
    spectrumZoomPanel.add(mz_zoomAll_,new GridBagConstraints(0, 3, 2, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(3, 0, 0, 0), 0, 0));
    
    exportSpectra_ = new ExportPanel(null,Color.BLACK,this,false,false,true);
    spectrumZoomPanel.add(exportSpectra_,new GridBagConstraints(0, 4, 2, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(3, 0, 0, 0), 0, 0));
    
    if (exportChromatogramsFromDRView_)
    {
    	l2dMenu.add(exportSpectra_,new GridBagConstraints(0, 11, 2, 1, 0.0, 0.0
          ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(15, 0, 0, 0), 0, 0));
    }
    
    exportFileChooser_ = new JFileChooser();
    exportFileChooser_.setPreferredSize(new Dimension(600,500));
    
  }
  
  private void initHelpPanel(){
    helpPanel_.setLayout(new BorderLayout());
    JPanel centerPanel = new JPanel();
    helpPanel_.add(centerPanel,BorderLayout.CENTER);
    centerPanel.setLayout(new GridBagLayout());
    
    JLabel helpText = new JLabel("To access the user manual on click one of the links:");
    helpText.setToolTipText(TooltipTexts.HELP_USER_MANUAL);
    centerPanel.add(helpText,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    String releaseWOMinor = Settings.VERSION.substring(0,Settings.VERSION.lastIndexOf("."));
    String linkAddress = "http://genome.tugraz.at/lda2/"+releaseWOMinor+"/LDA_"+releaseWOMinor+".pdf";

    
    JHyperlink linkText = new JHyperlink(linkAddress,linkAddress);
    centerPanel.add(linkText,new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    
    String baseRootPath = getClass().getResource("/at/tugraz/genome/lda/LipidDataAnalyzer.class").getPath();
    baseRootPath = baseRootPath.substring(0,baseRootPath.length()-"at/tugraz/genome/lda/LipidDataAnalyzer.class".length());
    if (baseRootPath.indexOf("!")>-1) baseRootPath = baseRootPath.substring(0,baseRootPath.lastIndexOf("!"));
    baseRootPath = baseRootPath.substring(0,baseRootPath.lastIndexOf("/"));
    String rootPath = baseRootPath+"/doc/LDA.pdf";
//      rootPath = (new File((new URL(rootPath)).toURI())).getCanonicalPath();
    rootPath = rootPath.replaceAll("\\\\", "/");
      
    if (rootPath.indexOf("file:")!=-1)rootPath = rootPath.substring(rootPath.indexOf("file:")+"file:".length());
    while (rootPath.startsWith("/")) rootPath = rootPath.substring(1);    
    rootPath = "file:///"+rootPath;
    linkText = new JHyperlink("Local instance in /doc folder",rootPath);
    centerPanel.add(linkText,new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    centerPanel.add(new JLabel(" "),new GridBagConstraints(0, 3, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));

    
    helpText = new JLabel("To access the sphingolipid data description document click one of the links:");
    helpText.setToolTipText(TooltipTexts.HELP_EXAMPLES);
    centerPanel.add(helpText,new GridBagConstraints(0, 4, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    linkAddress = "http://genome.tugraz.at/lda2/data/DataDescriptionSphingolipids.pdf";
    linkText = new JHyperlink(linkAddress,linkAddress);
    centerPanel.add(linkText,new GridBagConstraints(0, 5, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    
    rootPath = baseRootPath+"/examples/SphingolipidDataDescription.pdf";
//      rootPath = (new File((new URL(rootPath)).toURI())).getCanonicalPath();
    rootPath = rootPath.replaceAll("\\\\", "/");
    if (rootPath.indexOf("file:")!=-1)rootPath = rootPath.substring(rootPath.indexOf("file:")+"file:".length());
    while (rootPath.startsWith("/")) rootPath = rootPath.substring(1);    
    rootPath = "file:///"+rootPath;
    linkText = new JHyperlink("Local instance in /examples folder",rootPath);
    centerPanel.add(linkText,new GridBagConstraints(0, 6, 1, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    
    centerPanel.add(new JLabel(" "),new GridBagConstraints(0, 7, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));

    
    helpText = new JLabel("To access the phospholipid data description document click one of the links:");
    helpText.setToolTipText(TooltipTexts.HELP_EXAMPLES);
    centerPanel.add(helpText,new GridBagConstraints(0, 8, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    linkAddress = "http://genome.tugraz.at/lda2/data/DataDescription.pdf";
    linkText = new JHyperlink(linkAddress,linkAddress);
    centerPanel.add(linkText,new GridBagConstraints(0, 9, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    
    rootPath = baseRootPath+"/examples/DataDescription.pdf";
//      rootPath = (new File((new URL(rootPath)).toURI())).getCanonicalPath();
    rootPath = rootPath.replaceAll("\\\\", "/");
      
    if (rootPath.indexOf("file:")!=-1)rootPath = rootPath.substring(rootPath.indexOf("file:")+"file:".length());
    while (rootPath.startsWith("/")) rootPath = rootPath.substring(1);    
    rootPath = "file:///"+rootPath;
    linkText = new JHyperlink("Local instance in /examples folder",rootPath);
    centerPanel.add(linkText,new GridBagConstraints(0, 10, 1, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    
    centerPanel.add(new JLabel(" "),new GridBagConstraints(0, 11, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    helpText = new JLabel("The study data is available from our homepage:");
    helpText.setToolTipText(TooltipTexts.HELP_EXAMPLE_DOWNLOAD);
    centerPanel.add(helpText,new GridBagConstraints(0, 12, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    linkText = new JHyperlink("Link to Genome download","http://genome.tugraz.at/lda2/lda_data.shtml");
    centerPanel.add(linkText,new GridBagConstraints(0, 13, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    helpText = new JLabel("Please consult the examples document chapter to avoid unnecessary downloads!");
    centerPanel.add(helpText,new GridBagConstraints(0, 14, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));

    centerPanel.add(new JLabel(" "),new GridBagConstraints(0, 15, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));    
  }
  
  //TODO: Add uni graz logo? put this in it's own singleton class? Add Uni Graz website?
  private void initAboutPanel(){
    aboutPanel_.setLayout(new BorderLayout());
    JPanel topPanel = new JPanel();
    aboutPanel_.add(topPanel,BorderLayout.NORTH);
    topPanel.setLayout(new BorderLayout());
    JPanel logoPanel = new JPanel();
    logoPanel.setLayout(new GridBagLayout());
    topPanel.add(logoPanel,BorderLayout.WEST);
    JPanel headerPanel = new JPanel();
//    topPanel.add(headerPanel,BorderLayout.CENTER);
    logoPanel.add(new JLabel(ldaLogo_), new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(20, 50, 0, 0), 0, 0));
    logoPanel.add(new JLabel(tugLogo_), new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(10, 50, 0, 0), 0, 0));
    logoPanel.add(headerPanel, new GridBagConstraints(1, 0, 1, 2, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(10, 50, 0, 0), 0, 0));

    
//    logoPanel.setLayout(new BorderLayout());
//    logoPanel.add(new JLabel(ldaLogo_),BorderLayout.CENTER);
//    JLabel tugLabel = new JLabel(tugLogo_);
//    logoPanel.add(tugLabel,BorderLayout.SOUTH);
    
    
    headerPanel.setLayout(new GridBagLayout());
    JLabel text = new JLabel("L D A");
    text.setFont(new Font("Arial",Font.BOLD, 50));
    headerPanel.add(text,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    text = new JLabel("Lipid Data Analyzer");
    text.setFont(new Font("Arial",Font.BOLD, 24));
    headerPanel.add(text,new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    headerPanel.add(new JLabel(" "),new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    text = new JLabel("Designed and developed by:");
    text.setFont(new Font("Arial",Font.PLAIN, 18));
    headerPanel.add(text,new GridBagConstraints(0, 3, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    text = new JLabel("J\u00fcrgen Hartler, Alexander Triebl, Martin Tr\u00f6tzm\u00fcller, Andreas Ziegl, Leonida M. Lamp");
    text.setFont(new Font("Arial",Font.PLAIN, 16));
    headerPanel.add(text,new GridBagConstraints(0, 4, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    headerPanel.add(new JLabel(" "),new GridBagConstraints(0, 5, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    text = new JLabel("Thallinger Lab");
    text.setFont(new Font("Arial",Font.PLAIN, 16));
    headerPanel.add(text,new GridBagConstraints(0, 6, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    //Removed according to the instruction of the corresponding author: Gerhard G Thallinger
//    text = new JLabel("Institute of Computational Biotechnology");
//    text.setFont(new Font("Arial",Font.PLAIN, 16));
//    headerPanel.add(text,new GridBagConstraints(0, 7, 1, 1, 0.0, 0.0
//        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    text = new JLabel("Graz University of Technology");
    text.setFont(new Font("Arial",Font.PLAIN, 16));
    headerPanel.add(text,new GridBagConstraints(0, 8, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    JHyperlink linkText = new JHyperlink("http://genome.tugraz.at/lda2","http://genome.tugraz.at/lda2");
    linkText.setFont(new Font("Arial",Font.PLAIN, 16));
    headerPanel.add(linkText,new GridBagConstraints(0, 9, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    
    text = new JLabel("Lipid Data Analyzer");
    text.setFont(new Font("Arial",Font.PLAIN, 12));
    logoPanel.add(text,new GridBagConstraints(0, 3, 3, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(30, 50, 0, 0), 0, 0));
    text = new JLabel("Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry");
    text.setFont(new Font("Arial",Font.PLAIN, 12));
    logoPanel.add(text,new GridBagConstraints(0, 4, 3, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 50, 0, 0), 0, 0));
    text = new JLabel("Copyright \u00A9 2024 J\u00fcrgen Hartler, Andreas Ziegl, Gerhard G Thallinger, Leonida M Lamp");
    text.setFont(new Font("Arial",Font.PLAIN, 12));
    logoPanel.add(text,new GridBagConstraints(0, 5, 3, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 50, 0, 0), 0, 0));
    text = new JLabel("This program is free software: you can redistribute it and/or modify");
    text.setFont(new Font("Arial",Font.PLAIN, 12));
    logoPanel.add(text,new GridBagConstraints(0, 6, 2, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(30, 50, 0, 0), 0, 0));
    text = new JLabel("it under the terms of the GNU General Public License as published by");
    text.setFont(new Font("Arial",Font.PLAIN, 12));
    logoPanel.add(text,new GridBagConstraints(0, 7, 2, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 50, 0, 0), 0, 0));
    text = new JLabel("the Free Software Foundation, either version 3 of the License, or");
    text.setFont(new Font("Arial",Font.PLAIN, 12));
    logoPanel.add(text,new GridBagConstraints(0, 8, 2, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 50, 0, 0), 0, 0));
    text = new JLabel("(at your option) any later version.");
    text.setFont(new Font("Arial",Font.PLAIN, 12));
    logoPanel.add(text,new GridBagConstraints(0, 9, 2, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 50, 0, 0), 0, 0));

    text = new JLabel("This program is distributed in the hope that it will be useful,");
    text.setFont(new Font("Arial",Font.PLAIN, 12));
    logoPanel.add(text,new GridBagConstraints(0, 10, 2, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(30, 50, 0, 0), 0, 0));
    text = new JLabel("but WITHOUT ANY WARRANTY; without even the implied warranty of");
    text.setFont(new Font("Arial",Font.PLAIN, 12));
    logoPanel.add(text,new GridBagConstraints(0, 11, 2, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 50, 0, 0), 0, 0));
    text = new JLabel("MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the");
    text.setFont(new Font("Arial",Font.PLAIN, 12));
    logoPanel.add(text,new GridBagConstraints(0, 12, 2, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 50, 0, 0), 0, 0));
    text = new JLabel("GNU General Public License for more details.");
    text.setFont(new Font("Arial",Font.PLAIN, 12));
    logoPanel.add(text,new GridBagConstraints(0, 13, 2, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 50, 0, 0), 0, 0));

    text = new JLabel(" You should have received a copy of the GNU General Public License");
    text.setFont(new Font("Arial",Font.PLAIN, 12));
    logoPanel.add(text,new GridBagConstraints(0, 14, 2, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(30, 50, 0, 0), 0, 0));
    JPanel linkInText = new JPanel();
    linkInText.setLayout(new GridBagLayout());

    text = new JLabel("along with this program.  If not, see <");
    text.setFont(new Font("Arial",Font.PLAIN, 12));
    linkInText.add(text,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    linkInText.add(text);
    linkText = new JHyperlink("http://www.gnu.org/licenses/","http://www.gnu.org/licenses/");
    linkText.setFont(new Font("Arial",Font.PLAIN, 12));
    linkInText.add(linkText,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    text = new JLabel(">.");
    text.setFont(new Font("Arial",Font.PLAIN, 12));
    linkInText.add(text,new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    
    logoPanel.add(linkInText,new GridBagConstraints(0, 15, 2, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 50, 0, 0), 0, 0));
    

    if (Settings.isWindows()){
      text = new JLabel("For msconvert file translation (ProteoWizard), the following third-party software is provided:");
      text.setFont(new Font("Arial",Font.PLAIN, 12));
      logoPanel.add(text,new GridBagConstraints(0, 16, 2, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(30, 50, 0, 0), 0, 0));

      text = new JLabel("WIFF Reader Distributable Beta SDK. Copyright \u00A9 2013 AB SCIEX");
      text.setFont(new Font("Arial",Font.PLAIN, 12));
      logoPanel.add(text,new GridBagConstraints(0, 17, 2, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(5, 50, 0, 0), 0, 0));
      text = new JLabel("MASSHUNTER DATA ACCESS COMPONENT RUNTIME VERSION. Copyright \u00A9 2016 Agilent Technologies");
      text.setFont(new Font("Arial",Font.PLAIN, 12));
      logoPanel.add(text,new GridBagConstraints(0, 18, 2, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 50, 0, 0), 0, 0));
      text = new JLabel("This software uses CompassXtract software. Copyright \u00A9 2011, 2013, 2013 by Bruker Daltonik GmbH. All rights reserved.");
      text.setFont(new Font("Arial",Font.PLAIN, 12));
      logoPanel.add(text,new GridBagConstraints(0, 19, 2, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 50, 0, 0), 0, 0));
      text = new JLabel("MSFileReader file reading tool. Copyright \u00A9 2009 - 2014 by Thermo Fisher Scientific, Inc. All rights reserved.");
      text.setFont(new Font("Arial",Font.PLAIN, 12));
      logoPanel.add(text,new GridBagConstraints(0, 20, 2, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 50, 0, 0), 0, 0));
      text = new JLabel("These are individual software packages which in no event represent combined work with Lipid Data Analyzer.");
      text.setFont(new Font("Arial",Font.PLAIN, 12));
      logoPanel.add(text,new GridBagConstraints(0, 21, 2, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(5, 50, 0, 0), 0, 0));
      text = new JLabel("None of these software are governed by the license the Lipid Data Analyzer (GNU GPL) is providing you.");
      text.setFont(new Font("Arial",Font.PLAIN, 12));
      logoPanel.add(text,new GridBagConstraints(0, 22, 2, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 50, 0, 0), 0, 0));
    }

  }

  private void initSettingsPanel(){
    settingsPanel_.setLayout(new BorderLayout());
    JPanel centerPanel = new JPanel();
    settingsPanel_.add(centerPanel,BorderLayout.CENTER);
    centerPanel.setLayout(new GridBagLayout());
    centerPanel.add(new JLabel("Please do not change the settings while a calculation is running!"), new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0,0,0,0), 0, 0));

    
    JPanel inputPanel = new JPanel();
    inputPanel.setLayout(new GridBagLayout());
    JLabel msMachineLabel = new JLabel("MS settings: ");
    msMachineLabel.setToolTipText(TooltipTexts.SETTINGS_MS_MACHINE);
    inputPanel.add(msMachineLabel, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(1, 1, 1, 1), 0, 0));
    msMachineTypes_ = new JComboBox<String>();
    for (String msSetting : Settings.getPropertyFileNames()) msMachineTypes_.addItem(msSetting);
    String currentMachine = LipidomicsConstants.getCurrentMSMachine();
    msMachineTypes_.setSelectedItem(currentMachine);
    msMachineTypes_.setToolTipText(TooltipTexts.SETTINGS_MS_MACHINE);
    msMachineTypes_.addItemListener(new FragSettingsChangeListener("ChangeMachine"));
    inputPanel.add(msMachineTypes_, new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(1, 1, 1, 1), 0, 0));
    centerPanel.add(inputPanel, new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0,0,0,0), 0, 0));
        
    centerPanel.add(new JLabel(" "), new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0,0,0,0), 0, 0));
    centerPanel.add(new JLabel("Please do not change the fragmentation settings while a calculation is running!"), new GridBagConstraints(0, 3, 1, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0,0,0,0), 0, 0));

    JPanel inputFragSettings1 = new JPanel();
    inputFragSettings1.setLayout(new GridBagLayout());
    JLabel fragmentaionLabel1 = new JLabel("Fragmentation Selection 1: ");
    fragmentaionLabel1.setToolTipText(TooltipTexts.SETTINGS_MSN_FRAGMENTATION);
    inputFragSettings1.add(fragmentaionLabel1, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(1, 1, 1, 1), 0, 0));
    fragmentationSettings1_ = new JComboBox<String>();
    fragmentationSettings1_.setToolTipText(TooltipTexts.SETTINGS_MSN_FRAGMENTATION);
    inputFragSettings1.add(fragmentationSettings1_, new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(1, 1, 1, 1), 0, 0));
    
    inputFragSettings2_ = new JPanel();
    inputFragSettings2_.setLayout(new GridBagLayout());
    JLabel fragmentaionLabel2 = new JLabel("Fragmentation Selection 2: ");
    fragmentaionLabel2.setToolTipText(TooltipTexts.SETTINGS_MSN_FRAGMENTATION);
    inputFragSettings2_.add(fragmentaionLabel2, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(1, 1, 1, 1), 0, 0));
    fragmentationSettings2_ = new JComboBox<String>();
    fragmentationSettings2_.setToolTipText(TooltipTexts.SETTINGS_MSN_FRAGMENTATION);
    inputFragSettings2_.add(fragmentationSettings2_, new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(1, 1, 1, 1), 0, 0));
    
    refreshFragSettingsSelection();
    fragmentationSettings1_.addItemListener(new FragSettingsChangeListener("ChangeFragSelection"));
    centerPanel.add(inputFragSettings1, new GridBagConstraints(0, 4, 1, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0,0,0,0), 0, 0));
    centerPanel.add(inputFragSettings2_, new GridBagConstraints(0, 5, 1, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0,0,0,0), 0, 0));

    JPanel buttonPanel = new JPanel();
    JButton applyButton = new JButton("Apply");
    applyButton.setActionCommand("ApplyOtherMachineSettings");
    applyButton.setToolTipText(TooltipTexts.SETTINGS_BUTTON_APPLY);
    applyButton.addActionListener(this);
    buttonPanel.add(applyButton);
    JButton saveButton = new JButton("Save as default");
    saveButton.setActionCommand("SaveOtherMachineSettings");
    saveButton.setToolTipText(TooltipTexts.SETTINGS_BUTTON_SAVE);
    saveButton.addActionListener(this);
    buttonPanel.add(saveButton);
    centerPanel.add(buttonPanel, new GridBagConstraints(0, 6, 1, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));

  }
    
  private void initTimer(){
    timer_ = new java.util.Timer();
    timer_.schedule(new ThreadSupervisor(), 10, 1000);
  }
  
  private void createDisplayTopMenu(){
    this.displayTopMenu = new JPanel();
    this.displayTopMenu.setLayout(new GridBagLayout());
    
    this.selectedChromFile = new JTextField(62);
    selectedChromFile.setMinimumSize(selectedChromFile.getPreferredSize());
    selectedChromFile.setToolTipText(TooltipTexts.DISPLAY_OPEN_CHROM);
    displayTopMenu.add(selectedChromFile,new GridBagConstraints(0, 0, 6, 1, 0.0, 0.0
      ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    jButtonChromOpen = new JButton("Open Chrom");
    jButtonChromOpen.addActionListener(this);
    jButtonChromOpen.setActionCommand("showChromFileChooser");
    jButtonChromOpen.setToolTipText(TooltipTexts.DISPLAY_OPEN_CHROM);
    displayTopMenu.add(jButtonChromOpen,new GridBagConstraints(7, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 6, 0, 0), 0, 0));
    this.selectedResultFile = new JTextField(62);
    selectedResultFile.setMinimumSize(selectedResultFile.getPreferredSize());
    selectedResultFile.setToolTipText(TooltipTexts.DISPLAY_OPEN_RESULT);
    displayTopMenu.add(selectedResultFile,new GridBagConstraints(0, 1, 6, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    this.jButtonResultOpen = new JButton("Open Result");
    jButtonResultOpen.addActionListener(this);
    jButtonResultOpen.setActionCommand("showResultChooser");
    jButtonResultOpen.setToolTipText(TooltipTexts.DISPLAY_OPEN_RESULT);
    displayTopMenu.add(jButtonResultOpen,new GridBagConstraints(7, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 6, 0, 0), 0, 0));
    startDisplay = new JButton("Start Display");
    startDisplay.addActionListener(this);
    startDisplay.setActionCommand("startDisplay");
    startDisplay.setToolTipText(TooltipTexts.DISPLAY_START);
    displayTopMenu.add(startDisplay,new GridBagConstraints(8, 0, 1, 2, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0)); 
    
    jButtonResultOpen.setPreferredSize(jButtonChromOpen.getPreferredSize());
  }
  
  private void createResultsMenu(){
    resultsMenu_ = new JPanel();
    resultsMenu_.setLayout(new GridBagLayout());
    
    
    JPanel buttonPanel = new JPanel();
    buttonPanel.setLayout(new GridBagLayout()); 
    
    resultsMenu_.add(buttonPanel,new GridBagConstraints(0, 0, 5, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));    
    
    jButtonResultFilesOpen_ = new JButton("Add Files", addFilesIcon_);
    jButtonResultFilesOpen_.addActionListener(this);
    jButtonResultFilesOpen_.setActionCommand("showResultsFilesChooser");
    jButtonResultFilesOpen_.setToolTipText(TooltipTexts.STATISTICS_ADD_RESULT_FILES);
    buttonPanel.add(jButtonResultFilesOpen_,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 5, 2, 0), 0, 0));
    
    jButtonResultsDirOpen_ = new JButton("Add Dir",addFolderIcon_);
    jButtonResultsDirOpen_.addActionListener(this);
    jButtonResultsDirOpen_.setActionCommand("showResultsDirChooser");
    jButtonResultsDirOpen_.setToolTipText(TooltipTexts.STATISTICS_ADD_DIR);
    buttonPanel.add(jButtonResultsDirOpen_,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 2, 2, 0), 0, 0));
    
    jButtonResultFilesRemove_ = new JButton("Remove",removeFilesIcon_);
    jButtonResultFilesRemove_.addActionListener(this);
    jButtonResultFilesRemove_.setActionCommand("removeResultFiles");
    jButtonResultFilesRemove_.setToolTipText(TooltipTexts.STATISTICS_REMOVE_SELECTION);
    buttonPanel.add(jButtonResultFilesRemove_,new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 2, 2, 0), 0, 0));    
    
    jButtonResultFilesClean_ = new JButton("Remove all",removeFilesIcon_);
    jButtonResultFilesClean_.addActionListener(this);
    jButtonResultFilesClean_.setActionCommand("removeAllResultFiles");
    jButtonResultFilesClean_.setToolTipText(TooltipTexts.STATISTICS_REMOVE_ALL_SELECTION);
    buttonPanel.add(jButtonResultFilesClean_,new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 2, 2, 0), 0, 0));
    
    jButtonResultAddToGroup_ = new JButton("Add to group",addToGroupIcon_);
    jButtonResultAddToGroup_.addActionListener(this);
    jButtonResultAddToGroup_.setActionCommand("addToGroup");
    jButtonResultAddToGroup_.setToolTipText(TooltipTexts.STATISTICS_ADD_TO_GROUP);
    buttonPanel.add(jButtonResultAddToGroup_,new GridBagConstraints(4, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 2, 2, 0), 0, 0));    

    analysisSelectionTablePanel_ = new JPanel();
    this.generateResultsAnalysisTablePane(new String[0][0]);
//    resultsMenu_.add(analysisSelectionTablePanel_);
    resultsMenu_.add(analysisSelectionTablePanel_,new GridBagConstraints(0, 1, 5, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));  
    
    JPanel statisticalAnalysisSettingsPanel = new JPanel();
    statisticalAnalysisSettingsPanel.setLayout(new GridBagLayout()); 
    
    //frameHeight_+5 is the resulting dimension when analysisTablePane is added to analysisSelectionTablePanel_
    statisticalAnalysisSettingsPanel.setPreferredSize(new Dimension(FRAME_HEIGHT+5, 350));
    TitledBorder title;
    Border loweredEtched = BorderFactory.createEtchedBorder(EtchedBorder.LOWERED);
    title = BorderFactory.createTitledBorder(loweredEtched,"Settings for Statistical Analysis");
    statisticalAnalysisSettingsPanel.setBorder(title);
    
    groupsPanel_= new GroupsPanel();
    resultsMenu_.add(groupsPanel_ ,new GridBagConstraints(0, 2, 5, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    
    cutoffSettingsPanel_ = new CutoffSettingsPanel(TooltipTexts.STATISTICS_ADD_CUTOFF_SETTINGS);
    resultsMenu_.add(cutoffSettingsPanel_ ,new GridBagConstraints(0, 3, 5, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));  
    
    quantSettingsPanel_ = new AbsoluteQuantSettingsPanel (groupsPanel_);
    resultsMenu_.add(quantSettingsPanel_ ,new GridBagConstraints(0, 4, 5, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0)); 
    
    resultsMenu_.add(statisticalAnalysisSettingsPanel,new GridBagConstraints(0, 5, 5, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(30, 0, 0, 0), 0, 0));
    
    JPanel groupRtPanel = new JPanel();
    groupRtPanel.setLayout(new GridBagLayout());  
    JLabel sepRtLabel = new JLabel("Show hits with different RT separately: ");
    sepRtLabel.setToolTipText(TooltipTexts.STATISTICS_SEPARATE_RT);
    groupRtPanel.add(sepRtLabel,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    
    separateHitsByRT_  = new JCheckBox();
    separateHitsByRT_.setSelected(LipidomicsConstants.isShotgun()!=1);
    separateHitsByRT_.setActionCommand(CHANGE_SEPARATE_RT_STATUS);
    separateHitsByRT_.addActionListener(this);
    separateHitsByRT_.setToolTipText(TooltipTexts.STATISTICS_SEPARATE_RT);
    groupRtPanel.add(separateHitsByRT_,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    
    rtGroupingTime_ = new JTextField(3);
    rtGroupingTime_.setMinimumSize(rtGroupingTime_.getPreferredSize());
    rtGroupingTime_.setHorizontalAlignment(JTextField.RIGHT);
    rtGroupingTime_.setInputVerifier(new DoubleVerifier(true));
    rtGroupingTime_.setText("0.1");
    rtGroupingTime_.setToolTipText(TooltipTexts.STATISTICS_SEPARATE_RT);
    groupRtPanel.add(rtGroupingTime_,new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    rtGroupingTime_.setEnabled(separateHitsByRT_.isSelected());
    
    rtTimeUnit_ = new JLabel("min");
    rtTimeUnit_.setToolTipText(TooltipTexts.STATISTICS_SEPARATE_RT);
    groupRtPanel.add(rtTimeUnit_,new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    rtTimeUnit_.setEnabled(separateHitsByRT_.isSelected());
    
    JSeparator js = new JSeparator(SwingConstants.VERTICAL);
    js.setPreferredSize(new Dimension(5,1));
    GridBagConstraints gbc = new GridBagConstraints(4, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0);
    gbc.fill = GridBagConstraints.VERTICAL;
    gbc.weighty = 1;
    groupRtPanel.add(js, gbc);
    
    showMSnEvidenceStat_ = new JCheckBox("Load MSn Only");
    showMSnEvidenceStat_.setToolTipText(TooltipTexts.STATISTICS_SHOW_MSN_ONLY);
    showMSnEvidenceStat_.addItemListener(new LipidomicsItemListener("showMSnOnlyStat"));
    groupRtPanel.add(showMSnEvidenceStat_ ,new GridBagConstraints(5, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    
    
    js = new JSeparator(SwingConstants.VERTICAL);
    js.setPreferredSize(new Dimension(5,1));
    gbc = new GridBagConstraints(6, 0, 1, 1, 0.0, 1.0
            ,GridBagConstraints.WEST, GridBagConstraints.VERTICAL, new Insets(0, 6, 0, 0), 0, 0);
    groupRtPanel.add(js, gbc);
    
    showChainEvidenceStat_ = new JCheckBox("Chain Evidence Only");
    showChainEvidenceStat_.setEnabled(false);
    showChainEvidenceStat_.setToolTipText(TooltipTexts.STATISTICS_CHAIN_EVIDENCE_ONLY);
    showChainEvidenceStat_.addItemListener(new LipidomicsItemListener("showChainOnlyStat"));
    groupRtPanel.add(showChainEvidenceStat_ ,new GridBagConstraints(7, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    
    js = new JSeparator(SwingConstants.VERTICAL);
    js.setPreferredSize(new Dimension(5,1));
    gbc.gridx = 8;
    groupRtPanel.add(js, gbc);
    
    combineWithOx_ = new JCheckBox("Combine Classes with ox-Classes");
    combineWithOx_.setToolTipText(TooltipTexts.STATISTICS_COMBINE_OX);
    combineWithOx_.addItemListener(new LipidomicsItemListener("combineWithOx"));
    groupRtPanel.add(combineWithOx_ ,new GridBagConstraints(9, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
 
    
    statisticalAnalysisSettingsPanel.add(groupRtPanel,new GridBagConstraints(0, 1, 5, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(5, 0, 15, 0), 0, 0));   
    
    JPanel correctOrderPanel = new JPanel();
    correctOrderPanel.setLayout(new GridBagLayout());  
    JLabel quantLabel = new JLabel("Quant file (for correct analyte order):");
    quantLabel.setToolTipText(TooltipTexts.STATISTICS_CORRECT_ORDER_FILE);
    correctOrderPanel.add(quantLabel,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    
    correctOrderFile_ = new JTextField(35);
    correctOrderFile_.setMinimumSize(correctOrderFile_.getPreferredSize());
    correctOrderFile_.setToolTipText(TooltipTexts.STATISTICS_CORRECT_ORDER_FILE);
    correctOrderPanel.add(correctOrderFile_,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
      ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    
    JButton correctOpen = new JButton("Select");
    correctOpen.addActionListener(this);
    correctOpen.setActionCommand("showCorrectOrderFileChooser");
    correctOpen.setToolTipText(TooltipTexts.STATISTICS_CORRECT_ORDER_FILE);
    correctOrderPanel.add(correctOpen,new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
      ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    
    if (Settings.useAlex()){
      JPanel correctOrderAlexPanel = new JPanel();
      correctOrderAlexPanel.setLayout(new GridBagLayout());
      correctOrderPanel.add(correctOrderAlexPanel,new GridBagConstraints(1, 1, 2, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
      JLabel ionModeLabel = new JLabel("Ion mode:");
      ionModeLabel.setToolTipText(TooltipTexts.QUANTITATION_ION_MODE);
      correctOrderPanel.add(ionModeLabel,new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
          ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
      ionModeOrder_ = new JComboBox<String>();
      ionModeOrder_.addItem("+");
      ionModeOrder_.addItem("-");
      ionModeOrder_.setToolTipText(TooltipTexts.QUANTITATION_ION_MODE);
      ionModeOrder_.setFont(SELECT_FIELD_FONT);
      correctOrderAlexPanel.add(ionModeOrder_,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
          ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
      JLabel ionModeInfo = new JLabel("(selection required for Alex123 target lists)");
      ionModeInfo.setToolTipText(TooltipTexts.QUANTITATION_ION_MODE);
      correctOrderAlexPanel.add(ionModeInfo,new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    }
    
    statisticalAnalysisSettingsPanel.add(correctOrderPanel,new GridBagConstraints(0, 2, 5, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 15, 0), 0, 0));    
    
    JLabel label = new JLabel("Internal-standard prefix: ");
    label.setToolTipText(TooltipTexts.STATISTICS_IS_PREFIX);
    statisticalAnalysisSettingsPanel.add(label ,new GridBagConstraints(0, 3, 2, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    
    internalStandardSelection_ = new JTextField(3);
    internalStandardSelection_.setMinimumSize(internalStandardSelection_.getPreferredSize());
    internalStandardSelection_.setText(Settings.getInternalStandardDefaultInput());
    internalStandardSelection_.setToolTipText(TooltipTexts.STATISTICS_IS_PREFIX);
    statisticalAnalysisSettingsPanel.add(internalStandardSelection_ ,new GridBagConstraints(2, 3, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    
    label = new JLabel("External-standard prefix: ");
    statisticalAnalysisSettingsPanel.add(label ,new GridBagConstraints(0, 4, 2, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    
    label.setToolTipText(TooltipTexts.STATISTICS_ES_PREFIX);
    externalStandardSelection_ = new JTextField(3);
    externalStandardSelection_.setMinimumSize(externalStandardSelection_.getPreferredSize());
    externalStandardSelection_.setText(Settings.getExternalStandardDefaultInput());
    externalStandardSelection_.setToolTipText(TooltipTexts.STATISTICS_ES_PREFIX);
    statisticalAnalysisSettingsPanel.add(externalStandardSelection_ ,new GridBagConstraints(2, 4, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));    
//    externalStandardSelection_.setToolTipText();
    
    JPanel cutoffPanel = new JPanel();
    cutoffLoadButton_ = new JButton("Load settings");
    cutoffLoadButton_.addActionListener(this);
    cutoffLoadButton_.setActionCommand("loadCutoffSettings");
    cutoffLoadButton_.setToolTipText(TooltipTexts.STATISTICS_CUTOFF_LOAD);
    cutoffPanel.add(cutoffLoadButton_);
    cutoffLoadButton_.setVisible(false);
    cutoffSaveButton_ = new JButton("Save settings");
    cutoffSaveButton_.addActionListener(this);
    cutoffSaveButton_.setActionCommand("saveCutoffSettings");
    cutoffSaveButton_.setToolTipText(TooltipTexts.STATISTICS_CUTOFF_SAVE);
    cutoffPanel.add(cutoffSaveButton_);
    cutoffSaveButton_.setVisible(false);
    saveCutoffSettingsFileChooser_ = new JFileChooser();
    saveCutoffSettingsFileChooser_.setPreferredSize(new Dimension(600,500));
    saveCutoffSettingsFileChooser_.setFileFilter(new FileNameExtensionFilter("Cutoff settings file (*.cut.xml)","xml"));
    loadCutoffSettingsFileChooser_ = new JFileChooser();
    loadCutoffSettingsFileChooser_.setPreferredSize(new Dimension(600,500));
    loadCutoffSettingsFileChooser_.setFileFilter(new FileNameExtensionFilter("Cutoff settings file (*.cut.xml)","xml"));
    jButtonResultCutoff_ = new JButton("Add cutoff settings");
    jButtonResultCutoff_.addActionListener(this);
    jButtonResultCutoff_.setActionCommand("addRemoveCutoffSettings");
    jButtonResultCutoff_.setToolTipText(TooltipTexts.STATISTICS_ADD_CUTOFF_SETTINGS);
    cutoffPanel.add(jButtonResultCutoff_);
    cutoffWarningLabel_ = new JLabel("(Please specify result-files before)");
    cutoffWarningLabel_.setToolTipText(TooltipTexts.STATISTICS_ADD_CUTOFF_SETTINGS);
    cutoffPanel.add(cutoffWarningLabel_); 
    statisticalAnalysisSettingsPanel.add(cutoffPanel,new GridBagConstraints(3, 3, 2, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));

    JPanel absQuantSetPanel = new JPanel();
    resultLoadButton_ = new JButton("Load settings");
    resultLoadButton_.addActionListener(this);
    resultLoadButton_.setActionCommand("loadAbsoluteSettings");
    resultLoadButton_.setToolTipText(TooltipTexts.STATISTICS_CUTOFF_LOAD);
    absQuantSetPanel.add(resultLoadButton_);
    resultLoadButton_.setVisible(false);
    resultSaveButton_ = new JButton("Save settings");
    resultSaveButton_.addActionListener(this);
    resultSaveButton_.setActionCommand("saveAbsoluteSettings");
    resultSaveButton_.setToolTipText(TooltipTexts.STATISTICS_CUTOFF_SAVE);
    absQuantSetPanel.add(resultSaveButton_);
    resultSaveButton_.setVisible(false);
    saveAbsSettingsFileChooser_ = new JFileChooser();
    saveAbsSettingsFileChooser_.setPreferredSize(new Dimension(600,500));
    saveAbsSettingsFileChooser_.setFileFilter(new FileNameExtensionFilter("Whole quant-settings file (*.wqs.xml)","xml"));
    loadAbsSettingsFileChooser_ = new JFileChooser();
    loadAbsSettingsFileChooser_.setPreferredSize(new Dimension(600,500));
    loadAbsSettingsFileChooser_.setFileFilter(new FileNameExtensionFilter("Whole quant-settings file (*.wqs.xml)","xml"));
    jButtonResultAbsQuant_ = new JButton("Add absolute settings");
    jButtonResultAbsQuant_.addActionListener(this);
    jButtonResultAbsQuant_.setActionCommand("addRemoveAbsoluteSettings");
    jButtonResultAbsQuant_.setToolTipText(TooltipTexts.STATISTICS_ADD_ABS_SETTINGS);
    absQuantSetPanel.add(jButtonResultAbsQuant_);
    resultWarningLabel_ = new JLabel("(Please specify result-files before)");
    resultWarningLabel_.setToolTipText(TooltipTexts.STATISTICS_ADD_ABS_SETTINGS);
    absQuantSetPanel.add(resultWarningLabel_); 
    statisticalAnalysisSettingsPanel.add(absQuantSetPanel,new GridBagConstraints(3, 4, 2, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    
    
    JButton maxSize = new JButton("Remove absolute settings"); //longest text size the following two buttons need to be adjusted to
    jButtonResultCutoff_.setPreferredSize(maxSize.getPreferredSize());
    jButtonResultAbsQuant_.setPreferredSize(maxSize.getPreferredSize());
    maxSize = null; //for garbage collection
    
    jButtonResultFilesAccept_ = new JButton("Accept");
    jButtonResultFilesAccept_.addActionListener(this);
    jButtonResultFilesAccept_.setActionCommand("acceptSelectedResultFiles");
    jButtonResultFilesAccept_.setToolTipText(TooltipTexts.STATISTICS_ACCEPT_SELECTION);
    statisticalAnalysisSettingsPanel.add(jButtonResultFilesAccept_,new GridBagConstraints(4, 7, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(40, 0, 0, 0), 0, 0));
    
//    JPanel omegaMassListPanel = new JPanel();
//    omegaMassListPanel.setLayout(new GridBagLayout()); 
//    
//    //x-Dimension of analysisTablePane + 5 as that's the resulting dimension when analysisTablePane is added to analysisSelectionTablePanel_
//    omegaMassListPanel.setPreferredSize(new Dimension(frameHeight_+5, 100));
//    title = BorderFactory.createTitledBorder(loweredEtched,"Create a MassList with retention times specific to \u03C9- double bond positions");
//    omegaMassListPanel.setBorder(title);
//    
//    resultsMenu_.add(omegaMassListPanel,new GridBagConstraints(0, 6, 5, 1, 0.0, 0.0
//        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
//    
//    jButtonOmegaMassList_ = new JButton("Create \u03C9-MassList");
//    jButtonOmegaMassList_.addActionListener(this);
//    jButtonOmegaMassList_.setActionCommand("createOmegaMassList");
//    jButtonOmegaMassList_.setToolTipText(TooltipTexts.STATISTICS_CREATE_MASSLIST);
//    omegaMassListPanel.add(jButtonOmegaMassList_,new GridBagConstraints(0, 0, 0, 0, 0.0, 0.0
//        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(10, 0, 0, 0), 0, 0));
    
  }
  
  public void actionPerformed(ActionEvent arg0)
  {
    String command = arg0.getActionCommand();
    if (command.equalsIgnoreCase("showMzxmlFileChooser")){
      if (this.mzXMLFileChooser==null)
        this.mzXMLFileChooser = new JFileChooser();
      this.mzXMLFileChooser.setPreferredSize(new Dimension(600,500));
      this.mzXMLFileChooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
      
      int returnVal = this.mzXMLFileChooser.showOpenDialog(new Frame());;
      if (returnVal == JFileChooser.APPROVE_OPTION) {
         String text = this.mzXMLFileChooser.getSelectedFile().getAbsolutePath();
         this.singleQuantMenu_.getSelectedMzxmlDirectory().setText(text);
      } else return;
    }
    if (command.equalsIgnoreCase("showMzxmlDirChooser")){
      if (this.mzXMLDirChooser==null)
        this.mzXMLDirChooser = new JFileChooser();
      this.mzXMLDirChooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
      this.mzXMLDirChooser.setPreferredSize(new Dimension(600,500));
      
      int returnVal = this.mzXMLDirChooser.showOpenDialog(new Frame());;
      if (returnVal == JFileChooser.APPROVE_OPTION) {
         String text = this.mzXMLDirChooser.getSelectedFile().getAbsolutePath();
         this.batchQuantMenu_.getSelectedMzxmlDirectory().setText(text);
      } else return;
    }
    String subCommand = "";
    if (command.equalsIgnoreCase("showResultsDirChooser")){
      if (this.resultsDirChooser_==null)
        this.resultsDirChooser_ = new JFileChooser();
      this.resultsDirChooser_.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
      this.resultsDirChooser_.setPreferredSize(new Dimension(600,500));
      
      int returnVal = this.resultsDirChooser_.showOpenDialog(new Frame());;
      if (returnVal == JFileChooser.APPROVE_OPTION) {
        Hashtable<String,File> avoidDuplicates = new Hashtable<String,File>();
        if (this.resultFiles_!=null){
          for (File file : this.resultFiles_) avoidDuplicates.put(file.getAbsolutePath(), file);
        }else{
          this.resultFiles_ = new Vector<File>();
        }
         String text = this.resultsDirChooser_.getSelectedFile().getAbsolutePath();
         File resultsDir = new File(text);
         if (resultsDir.exists() && resultsDir.isDirectory()){
           File[] resultFileCandidates = resultsDir.listFiles();
           for (int i=0; i!=resultFileCandidates.length;i++){
             if (resultFileCandidates[i].isFile() && !avoidDuplicates.containsKey(resultFileCandidates[i].getAbsolutePath())){
               String fileName = StaticUtils.extractFileName(resultFileCandidates[i].getAbsolutePath()); 
               String suffix = fileName.substring(fileName.lastIndexOf(".")+1);
               if ((suffix.equalsIgnoreCase("xls")||suffix.equalsIgnoreCase("xlsx")) && !fileName.startsWith(ExcelUtils.EXCEL_TEMP_PREFIX)){
                 avoidDuplicates.put(resultFileCandidates[i].getAbsolutePath(),resultFileCandidates[i]);
                 this.resultFiles_.add(resultFileCandidates[i]);
               }
             }
           }
           this.resultFiles_ = sortFilesByName(this.resultFiles_);
         }
         if (avoidDuplicates.size()>0){
           this.updateAnalysisSelectionTable();
         }
         if (jButtonResultAbsQuant_.getText().equalsIgnoreCase("Remove absolute settings")||jButtonResultCutoff_.getText().equalsIgnoreCase("Remove cutoff settings"))
           subCommand = "removeSettings";
         
      } else return;

    }
    if (command.equalsIgnoreCase("showResultsFilesChooser")){
      if (this.resultFilesChooser_==null){
        this.resultFilesChooser_ = new JFileChooser();
        this.resultFilesChooser_.setFileFilter(new FileNameExtensionFilter("LDA quant results (*.xls,*.xlsx)","xls","xlsx"));
      }  
      this.resultFilesChooser_.setMultiSelectionEnabled(true);
      this.resultFilesChooser_.setPreferredSize(new Dimension(600,500));  
      int returnVal = this.resultFilesChooser_.showOpenDialog(new Frame());;
      if (returnVal == JFileChooser.APPROVE_OPTION) {      
        File[] resultFileCandidates = this.resultFilesChooser_.getSelectedFiles();
        Hashtable<String,File> avoidDuplicates = new Hashtable<String,File>();
        if (this.resultFiles_!=null){
          for (File file : this.resultFiles_) avoidDuplicates.put(file.getAbsolutePath(), file);
        }else{
          this.resultFiles_ = new Vector<File>();
        }
        for (int i=0; i!=resultFileCandidates.length;i++){
          if (resultFileCandidates[i].isFile() && !avoidDuplicates.containsKey(resultFileCandidates[i].getAbsolutePath())){
            avoidDuplicates.put(resultFileCandidates[i].getAbsolutePath(),resultFileCandidates[i]);
            this.resultFiles_.add(resultFileCandidates[i]);
          }
        }
        if (avoidDuplicates.size()>0){
          this.updateAnalysisSelectionTable();
        }
        if (resultFileCandidates.length>0 && jButtonResultAbsQuant_.getText().equalsIgnoreCase("Remove absolute settings")||jButtonResultCutoff_.getText().equalsIgnoreCase("Remove cutoff settings"))
          subCommand = "removeSettings";
      } else return;
    }
    if (command.equalsIgnoreCase("showCorrectOrderFileChooser")){
      if (this.correctOrderFileChooser_==null){
        this.correctOrderFileChooser_ = new JFileChooser();
        if (Settings.useAlex()){
          this.correctOrderFileChooser_.setFileFilter(new FileNameExtensionFilter("LDA/ALEX quant file (*.xls,*.xlsx/*.txt)","xls","xlsx","txt"));
        } else
          this.correctOrderFileChooser_.setFileFilter(new FileNameExtensionFilter("LDA quant file (*.xls,*.xlsx)","xls","xlsx"));
      }  
      this.correctOrderFileChooser_.setPreferredSize(new Dimension(600,500));
      if (Settings.useAlex()) this.correctOrderFileChooser_.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
      else this.correctOrderFileChooser_.setFileSelectionMode(JFileChooser.FILES_ONLY);
      
      int returnVal = this.correctOrderFileChooser_.showOpenDialog(new Frame());;
      if (returnVal == JFileChooser.APPROVE_OPTION) {
         String text = this.correctOrderFileChooser_.getSelectedFile().getAbsolutePath();
         this.correctOrderFile_.setText(text);
      } else return;
    }
    if (this.resultFiles_!=null && command.equalsIgnoreCase("removeResultFiles")){
      int[] selectedColumns = resultFilesDisplayTable_.getSelectedRows();
      if (selectedColumns!=null && selectedColumns.length>0){
        List<Integer> selectedList = new ArrayList<Integer>();
        for (int i=0;i!=selectedColumns.length;i++){
          selectedList.add(selectedColumns[i]);
        }
        Collections.sort(selectedList);
        Vector<File> filesToRemove = new Vector<File>();
        for (int i=(selectedList.size()-1);i!=-1;i--){
          filesToRemove.add(this.resultFiles_.get((int)selectedList.get(i)));
          this.resultFiles_.remove((int)selectedList.get(i));
        }
        groupsPanel_.removeFiles(filesToRemove);
        this.updateAnalysisSelectionTable();
        if (jButtonResultAbsQuant_.getText().equalsIgnoreCase("Remove absolute settings")||jButtonResultCutoff_.getText().equalsIgnoreCase("Remove cutoff settings"))
          subCommand = "removeSettings";
      }
    }
    if (this.resultFiles_!=null && command.equalsIgnoreCase("removeAllResultFiles")){
      this.resultFiles_ = new Vector<File>();
      this.groupsPanel_.removeAllGroups();
      this.updateAnalysisSelectionTable();
      if (jButtonResultAbsQuant_.getText().equalsIgnoreCase("Remove absolute settings")||jButtonResultCutoff_.getText().equalsIgnoreCase("Remove cutoff settings"))
        subCommand = "removeSettings";
    }
    if (command.equalsIgnoreCase("addToGroup")){
      if (resultFiles_!=null){
        int[] selectedColumns = resultFilesDisplayTable_.getSelectedRows();
        if (selectedColumns!=null && selectedColumns.length>0){
          Vector<File> selectedFiles = new Vector<File>();
          for (int selectedColumn : selectedColumns){
            selectedFiles.add(resultFiles_.get(selectedColumn));
          }
          resultFilesDisplayTable_.getSelectionModel().clearSelection();
          int amountGroups = groupsPanel_.getGroups().size();
          InputDialog dlg = new InputDialog(new JFrame(), "Add to group", "Enter the group name", "Group "+String.valueOf(amountGroups+1));
          String groupName = dlg.getEnteredText();
          if (groupName!=null&&groupName.length()>0){
            groupName = groupName.trim();
            groupsPanel_.addGroup(groupName,selectedFiles,FRAME_HEIGHT);
          }else{
            new WarningMessage(new JFrame(), "Warning", "You have to specifiy a name for the group!");
          }
//        List<Integer> selectedList = new ArrayList<Integer>();
//        for (int i=0;i!=selectedColumns.length;i++){
//          selectedList.add(selectedColumns[i]);
//        }
//        Collections.sort(selectedList);
//        for (int i=(selectedList.size()-1);i!=-1;i--){
//          this.resultFiles_.remove((int)selectedList.get(i));
//        }
//        this.updateAnalysisSelectionTable();
        }else{
          new WarningMessage(new JFrame(), "Warning", "You have to select at least one result-file!");
        }
      }else{
        new WarningMessage(new JFrame(), "Warning", "First you have to add some result-files!");
      }
    }
    if (command.equalsIgnoreCase("acceptSelectedResultFiles")){
    	cleanupResultView();
//    	resultTabs_.setComponentAt(0, new LoadingPanel("Processing data, please wait..."));
//    	Thread thread = new Thread(new Runnable()
//  		{
//  			public void run()
//  			{
  				try
  				{
  					acceptResultFiles();
  				}
  				catch (ExcelInputFileException ex)
  				{
  					ex.printStackTrace();
  					new WarningMessage(new JFrame(), "Error", ex.getMessage());
  				}
//  				resultTabs_.setComponentAt(0, resultsSelectionPanel_);
//  				resultsSelectionPanel_.repaint(); //in case an error was thrown
//  			}
//  		});
//    	thread.start(); 
    }
    if (command.equalsIgnoreCase("showChromFileChooser")){
      if (chromFileChooser_==null)
        this.chromFileChooser_ = new JFileChooser();
      this.chromFileChooser_.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
      this.chromFileChooser_.setPreferredSize(new Dimension(600,500));
      
      int returnVal = this.chromFileChooser_.showOpenDialog(new Frame());
      if (returnVal == JFileChooser.APPROVE_OPTION) {
         String text = this.chromFileChooser_.getSelectedFile().getAbsolutePath();
         this.selectedChromFile.setText(text);

      } else return;
    }
    if (command.equalsIgnoreCase("showQuantFileChooser")){
      if (quantFileChooser_==null)
        this.quantFileChooser_ = new JFileChooser();
      this.quantFileChooser_.setPreferredSize(new Dimension(600,500));
      if (Settings.useAlex())
        this.quantFileChooser_.setFileFilter(new FileNameExtensionFilter("LDA/ALEX quant file (*.xls,*.xlsx/*.txt)","xls","xlsx","txt"));
      else this.quantFileChooser_.setFileFilter(new FileNameExtensionFilter("LDA quant file (*.xls,*.xlsx)","xls","xlsx"));
      
      int returnVal = this.quantFileChooser_.showOpenDialog(new Frame());
      if (returnVal == JFileChooser.APPROVE_OPTION) {
         String text = this.quantFileChooser_.getSelectedFile().getAbsolutePath();
         this.singleQuantMenu_.getSelectedQuantDir().setText(text);
      } else return;
    }
    if (command.equalsIgnoreCase("showQuantDirChooser")){
      if (quantDirChooser_==null)
        this.quantDirChooser_ = new JFileChooser();
      this.quantDirChooser_.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
      this.quantDirChooser_.setPreferredSize(new Dimension(600,500));
      
      int returnVal = this.quantDirChooser_.showOpenDialog(new Frame());;
      if (returnVal == JFileChooser.APPROVE_OPTION) {
         String text = this.quantDirChooser_.getSelectedFile().getAbsolutePath();
         this.batchQuantMenu_.getSelectedQuantDir().setText(text);
      } else return;
    }
    if (command.equalsIgnoreCase("showResultChooser")){
      if (resultFileChooser_==null){
        this.resultFileChooser_ = new JFileChooser();
        this.resultFileChooser_.setFileFilter(new FileNameExtensionFilter("LDA quant results (*.xls,*.xlsx)","xls","xlsx"));
      }  
      this.resultFileChooser_.setPreferredSize(new Dimension(600,500));
      int returnVal = this.resultFileChooser_.showOpenDialog(new Frame());;
      if (returnVal == JFileChooser.APPROVE_OPTION) {
         String text = this.resultFileChooser_.getSelectedFile().getAbsolutePath();
         this.selectedResultFile.setText(text);

      } else return;
    }
    if (command.equalsIgnoreCase("startBatchQuantification")){
      boolean canStartQuantification = false;
      float cutoff = 0f;
      float rtShift = 0f;
      if (batchQuantMenu_.getSelectedMzxmlDirectory().getText()!=null&&batchQuantMenu_.getSelectedMzxmlDirectory().getText().length()>0 &&
          batchQuantMenu_.getSelectedQuantDir().getText()!=null&&batchQuantMenu_.getSelectedQuantDir().getText().length()>0){
        canStartQuantification = true;
        if (this.batchQuantMenu_.getIsoValidation().isSelected()){
          if (this.batchQuantMenu_.getAmountOfIsotopes().getText()==null||this.batchQuantMenu_.getAmountOfIsotopes().getText().length()<1)
            canStartQuantification = false;
          else{
            try{
              Integer.parseInt(this.batchQuantMenu_.getAmountOfIsotopes().getText());
            }catch(NumberFormatException nfx){
              canStartQuantification = false;
            }
          }
          if (this.batchQuantMenu_.getAmountOfMatchingSearchIsotopes().getText()==null||this.batchQuantMenu_.getAmountOfMatchingSearchIsotopes().getText().length()<1)
            canStartQuantification = false;
          else{
            try{
              int isos = Integer.parseInt(this.batchQuantMenu_.getAmountOfMatchingSearchIsotopes().getText());
              if (isos<0)
                canStartQuantification = false;
            }catch(NumberFormatException nfx){
              canStartQuantification = false;
            }
          }
        }
        try{
          if (batchQuantMenu_.getCutoff().getText()!=null && batchQuantMenu_.getCutoff().getText().length()>0){
            cutoff = Float.parseFloat(batchQuantMenu_.getCutoff().getText().replaceAll(",", "."));
            LipidomicsConstants.getInstance().setRelativeMS1BasePeakCutoff(batchQuantMenu_.getCutoff().getText());
          }  
        }catch(NumberFormatException ex){new WarningMessage(new JFrame(), "Error", "The cutoff value must be float format!"); canStartQuantification=false;}
        try{
          if (batchQuantMenu_.getRtShift().getText()!=null && batchQuantMenu_.getRtShift().getText().length()>0)
            rtShift = Float.parseFloat(batchQuantMenu_.getRtShift().getText().replaceAll(",", "."));
        }catch(NumberFormatException ex){new WarningMessage(new JFrame(), "Error", "The RT-shift value must be float format!"); canStartQuantification=false;}

//        if (this.batchQuantMenu_.getSearchUnknownTime().isSelected()){
//        }
      }
      float minusTimeTol = 0f;
      float plusTimeTol = 0f;
      if (this.batchQuantMenu_.getTimeMinusTol().getText()!=null && this.batchQuantMenu_.getTimeMinusTol().getText().length()>0){
        try{ minusTimeTol = Float.parseFloat(this.batchQuantMenu_.getTimeMinusTol().getText());} catch (NumberFormatException nfx){new WarningMessage(new JFrame(), "Error", "The  \"Time before tol.\" value must be entered in float format!"); canStartQuantification=false;}
      }
      if (this.batchQuantMenu_.getTimePlusTol().getText()!=null && this.batchQuantMenu_.getTimePlusTol().getText().length()>0){
        try{ plusTimeTol = Float.parseFloat(this.batchQuantMenu_.getTimePlusTol().getText());} catch (NumberFormatException nfx){new WarningMessage(new JFrame(), "Error", "The  \"Time after tol.\" value must be entered in float format!"); canStartQuantification=false;}
      }
      if (canStartQuantification){      
        Vector<File> rawFiles = new Vector<File>();
        Vector<File> quantFiles = new Vector<File>();
        File rawDir = new File(this.batchQuantMenu_.getSelectedMzxmlDirectory().getText());
        File quantDir = new File(this.batchQuantMenu_.getSelectedQuantDir().getText());
        int amountOfIsotopes = 0;
        int isotopesMustMatch = 0;
        if (this.batchQuantMenu_.getIsoValidation().isSelected()){
          amountOfIsotopes = Integer.parseInt(this.batchQuantMenu_.getAmountOfIsotopes().getText());
          isotopesMustMatch = Integer.parseInt(this.batchQuantMenu_.getAmountOfMatchingSearchIsotopes().getText());
        }
//        if (this.batchQuantMenu_.getSearchUnknownTime().isSelected()){          
//        }
        if (rawDir.exists()&&rawDir.isDirectory()&&quantDir.exists()&&quantDir.isDirectory()){
          File[] rawFileCandidates = rawDir.listFiles();
          Hashtable<String,Vector<File>> avoidDuplication = new Hashtable<String,Vector<File>>();
          boolean mzXMLOrChromPresent = false;
          for (int i=0; i!=rawFileCandidates.length;i++){
            if (rawFileCandidates[i].isFile()){
              String[] fileNameAndSuffix = StaticUtils.extractFileNameAndSuffix(rawFileCandidates[i].getAbsolutePath()); 
              String suffix = fileNameAndSuffix[1];
              String fileName = fileNameAndSuffix[0];
              if (suffix.equalsIgnoreCase(AbstractXMLSpectraReader.FILE_TYPE_MZ_XML)||
              		suffix.equalsIgnoreCase(AbstractXMLSpectraReader.FILE_TYPE_MZ_ML)||
              		suffix.equalsIgnoreCase("raw")||
              		suffix.equalsIgnoreCase("chrom")||
              		suffix.equalsIgnoreCase("wiff")){
                if (suffix.equalsIgnoreCase(AbstractXMLSpectraReader.FILE_TYPE_MZ_XML)||suffix.equalsIgnoreCase(AbstractXMLSpectraReader.FILE_TYPE_MZ_ML)||suffix.equalsIgnoreCase("chrom")) mzXMLOrChromPresent = true;
                Vector<File> theFiles = new Vector<File>();
                if (avoidDuplication.containsKey(fileName)){
                  theFiles = avoidDuplication.get(fileName);
                }
                theFiles.add(rawFileCandidates[i]);
                avoidDuplication.put(fileName, theFiles);
              }
            }
            if (rawFileCandidates[i].isDirectory()){
              String[] fileNameAndSuffix = StaticUtils.extractFileNameAndSuffix(rawFileCandidates[i].getAbsolutePath()); 
              String suffix = fileNameAndSuffix[1];
              String fileName = fileNameAndSuffix[0];
              if (suffix.equalsIgnoreCase("raw")|| suffix.equalsIgnoreCase("d") ||suffix.equalsIgnoreCase("chrom")){
                if (suffix.equalsIgnoreCase("chrom")) mzXMLOrChromPresent = true;
                Vector<File> theFiles = new Vector<File>();
                if (avoidDuplication.containsKey(fileName)){
                  theFiles = avoidDuplication.get(fileName);
                }
                theFiles.add(rawFileCandidates[i]);
                avoidDuplication.put(fileName, theFiles);
              }
            }
          }
          for (String key : avoidDuplication.keySet()){
            Vector<File> theFiles = avoidDuplication.get(key);
            if (theFiles.size()==1){
              String suffix = StaticUtils.extractFileNameAndSuffix(theFiles.get(0).getAbsolutePath())[1];
              if (!mzXMLOrChromPresent || !suffix.equalsIgnoreCase("wiff"))
                rawFiles.add(theFiles.get(0));
            }else{
              int selectedIndex = -1;
              for (int i=0; i!=theFiles.size();i++){
                File file = theFiles.get(i);
                String suffix = file.getAbsolutePath().substring(file.getAbsolutePath().lastIndexOf(".")+1);
                if (mzXMLOrChromPresent && suffix.equalsIgnoreCase("wiff")) continue;
                if (suffix.equalsIgnoreCase("chrom")){
                  selectedIndex = i;
                }
              }
              if (selectedIndex>-1){
                rawFiles.add(theFiles.get(selectedIndex));
              }else{
                for (int i=0; i!=theFiles.size();i++){
                  File file = theFiles.get(i);
                  String suffix = file.getAbsolutePath().substring(file.getAbsolutePath().lastIndexOf(".")+1);
                  if (mzXMLOrChromPresent && suffix.equalsIgnoreCase("wiff")) continue;
                  if (suffix.equalsIgnoreCase(AbstractXMLSpectraReader.FILE_TYPE_MZ_XML) || suffix.equalsIgnoreCase(AbstractXMLSpectraReader.FILE_TYPE_MZ_ML)){
                    rawFiles.add(theFiles.get(i));
                  }
                }  
              }
            }
          }
          File[] quantificationFileCandidates = quantDir.listFiles();
          boolean containsTxtFiles = false;
          for (int i=0; i!=quantificationFileCandidates.length;i++){
            String suffix = quantificationFileCandidates[i].getAbsolutePath().substring(quantificationFileCandidates[i].getAbsolutePath().lastIndexOf(".")+1);
            if (suffix.equalsIgnoreCase("xls")||suffix.equalsIgnoreCase("xlsx")){
              quantFiles.add(quantificationFileCandidates[i]);
            } else if (suffix.equalsIgnoreCase("txt"))
              containsTxtFiles = true;
          }
          if (Settings.useAlex() && containsTxtFiles)
            quantFiles.add(quantDir);
          Vector<RawQuantificationPairVO> pairs = new Vector<RawQuantificationPairVO>();
          if (rawFiles.size()>0 && quantFiles.size()>0 && (pairs = generateQuantificationPairVOs(rawFiles,quantFiles)).size()>0){
            boolean ionMode = false;
            if (this.batchQuantMenu_.getIonMode()!=null && ((String)batchQuantMenu_.getIonMode().getSelectedItem()).equalsIgnoreCase("+"))
              ionMode = true;
            this.batchQuantMenu_.getBatchQuantTableModel().clearFiles();
            this.batchQuantMenu_.getBatchQuantTableModel().addFiles(pairs);
            this.batchQuantMenu_.getQuantifyingLabel().setText("Quantifying");
            this.batchQuantMenu_.getProgressBar().setValue(0);
            this.batchQuantMenu_.getQuantifyingPanel().setVisible(true);
            this.batchQuantMenu_.getStartQuantification().setEnabled(false);
            this.batchQuantMenu_.getSpinnerLabel().setVisible(true);
            batchQuantThread_ = new BatchQuantThread(batchQuantMenu_,
                minusTimeTol,plusTimeTol,amountOfIsotopes,isotopesMustMatch,cutoff, 
                rtShift,ionMode,false);
            batchQuantThread_.start();
          }else{
            if (rawFiles.size()==0){
              @SuppressWarnings("unused")
              WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "In the specified raw directory are no quantifyable files");
            }
            if (quantFiles.size()==0){
              @SuppressWarnings("unused")
              WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "In the specified quant directory are no quantifyable files");
            }
            if (rawFiles.size()>0 && quantFiles.size()>0)
              new WarningMessage(new JFrame(), "Warning", "In the specified directories are no quantifyable raw/quant pairs");
          }
        }else{
          if (!rawDir.exists()||!rawDir.isDirectory()){
            @SuppressWarnings("unused")
            WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "The raw directory does not exist");
          }
          if (!quantDir.exists()||!quantDir.isDirectory()){
            @SuppressWarnings("unused")
            WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "The quantification directory does not exist");
          }
        }
      }else{
//        boolean wasWarningMessage=false;
        if (batchQuantMenu_.getSelectedMzxmlDirectory().getText()==null||batchQuantMenu_.getSelectedMzxmlDirectory().getText().length()<1){
          @SuppressWarnings("unused")
          WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "You must specify a raw, mzXML, mzML or chrom file directory");
//          wasWarningMessage=true;
        }else if (batchQuantMenu_.getSelectedQuantDir().getText()==null||batchQuantMenu_.getSelectedQuantDir().getText().length()<1){
          @SuppressWarnings("unused")
          WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "You must specify a quant file directory");
//          wasWarningMessage=true;
        } else if (this.batchQuantMenu_.getIsoValidation().isSelected()){
          if (this.batchQuantMenu_.getAmountOfIsotopes().getText()==null||this.batchQuantMenu_.getAmountOfIsotopes().getText().length()<1){
            @SuppressWarnings("unused")
            WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "If you select the isotope option, please specify the amount of isotopes (integer format)!");
//            wasWarningMessage=true;
          }else{
            try{
              Integer.parseInt(this.batchQuantMenu_.getAmountOfIsotopes().getText());
            }catch(NumberFormatException nfx){
              @SuppressWarnings("unused")
              WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "The number of isotopes must be integer format!");
//              wasWarningMessage=true;
            }
          }
        }
        if (this.batchQuantMenu_.getAmountOfMatchingSearchIsotopes().getText()==null||this.batchQuantMenu_.getAmountOfMatchingSearchIsotopes().getText().length()<0){
          @SuppressWarnings("unused")
          WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "If you select the find molecules without retention-time option, please specify the amount of matching isotopes (integer format)!");
//          wasWarningMessage=true;
        }else{
          try{
            int isotopes = Integer.parseInt(this.batchQuantMenu_.getAmountOfMatchingSearchIsotopes().getText());
            if (isotopes<0){
              @SuppressWarnings("unused")
              WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "The number of matching isotopes must be >=0!");
//              wasWarningMessage=true;
            }
          }catch(NumberFormatException nfx){
            @SuppressWarnings("unused")
            WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "The number of matching isotopes must be integer format!");
//            wasWarningMessage=true;
          }
        } 

//        if (!wasWarningMessage && this.batchQuantMenu_.getSearchUnknownTime().isSelected()){
//        }
      }
    }
    if (command.equalsIgnoreCase("startQuantification")){
      boolean canStartQuantification = false;
      float beforeTolerance = 0f;
      float afterTolerance = 0f;
      float cutoff = 0f;
      float rtShift = 0f;
      if (singleQuantMenu_.getSelectedMzxmlDirectory().getText()!=null)
        singleQuantMenu_.getSelectedMzxmlDirectory().setText(singleQuantMenu_.getSelectedMzxmlDirectory().getText().trim());
      if (singleQuantMenu_.getSelectedQuantDir().getText()!=null)
        singleQuantMenu_.getSelectedQuantDir().setText(singleQuantMenu_.getSelectedQuantDir().getText().trim());
      if (singleQuantMenu_.getSelectedMzxmlDirectory().getText()!=null&&singleQuantMenu_.getSelectedMzxmlDirectory().getText().length()>0 &&
          singleQuantMenu_.getSelectedQuantDir().getText()!=null&&singleQuantMenu_.getSelectedQuantDir().getText().length()>0){
        canStartQuantification = true;
        if (!StaticUtils.existsFile(singleQuantMenu_.getSelectedMzxmlDirectory().getText())){
          new WarningMessage(new JFrame(), "Error", "The raw file \""+StaticUtils.extractFileName(singleQuantMenu_.getSelectedMzxmlDirectory().getText())+"\" does not exist!"); 
          return;
        }else{
          if (singleQuantMenu_.getSelectedQuantDir().getText().length()>3){
            String suffix = singleQuantMenu_.getSelectedQuantDir().getText().substring(singleQuantMenu_.getSelectedQuantDir().getText().lastIndexOf("."));
            if (!(suffix.equalsIgnoreCase(".xls")||(suffix.equalsIgnoreCase(".xlsx")||(Settings.useAlex()&&suffix.equalsIgnoreCase(".txt"))))){
              if (Settings.useAlex())
                new WarningMessage(new JFrame(), "Error", "For the mass lists just Excel files in the xls or xlsx or Text files in the txt format are allowed!"); 
              else
                new WarningMessage(new JFrame(), "Error", "For the mass lists just Excel files in the xls or xlsx format are allowed!"); 
              return;
            }
          }
        }
        if (canStartQuantification){
          String suffix = singleQuantMenu_.getSelectedMzxmlDirectory().getText().substring(singleQuantMenu_.getSelectedMzxmlDirectory().getText().lastIndexOf(".")+1);
          if (!(suffix.equalsIgnoreCase("chrom")||suffix.equalsIgnoreCase("head")||suffix.equalsIgnoreCase("idx")||
              suffix.equalsIgnoreCase("rtt")||suffix.equalsIgnoreCase("raw")||suffix.equalsIgnoreCase("d")||suffix.equalsIgnoreCase("wiff")||
              suffix.equalsIgnoreCase(AbstractXMLSpectraReader.FILE_TYPE_MZ_XML)||
              suffix.equalsIgnoreCase(AbstractXMLSpectraReader.FILE_TYPE_MZ_ML))){
            new WarningMessage(new JFrame(), "Error", "For the raw data just files and directories with the suffix raw, mzXML, mzML, d, wiff, and chrom are allowed!"); 
            return;
          }
        }
        if (!StaticUtils.existsFile(singleQuantMenu_.getSelectedQuantDir().getText())){
          new WarningMessage(new JFrame(), "Error", "The mass list file \""+StaticUtils.extractFileName(singleQuantMenu_.getSelectedQuantDir().getText())+"\" does not exist!"); 
          canStartQuantification=false;
        }
        try{
          if (singleQuantMenu_.getTimeMinusTol().getText()!=null && singleQuantMenu_.getTimeMinusTol().getText().length()>0)
            beforeTolerance = Float.parseFloat(singleQuantMenu_.getTimeMinusTol().getText().replaceAll(",", "."));
        }catch(NumberFormatException ex){new WarningMessage(new JFrame(), "Error", "The  \"Time before tol.\" value must be entered in float format!"); canStartQuantification=false;}
        try{
          if (singleQuantMenu_.getTimePlusTol().getText()!=null && singleQuantMenu_.getTimePlusTol().getText().length()>0)
            afterTolerance = Float.parseFloat(singleQuantMenu_.getTimePlusTol().getText().replaceAll(",", "."));
        }catch(NumberFormatException ex){new WarningMessage(new JFrame(), "Error", "The  \"Time after tol.\" value must be entered in float format!"); canStartQuantification=false;}
        if (this.singleQuantMenu_.getIsoValidation().isSelected()){
          if (this.singleQuantMenu_.getAmountOfIsotopes().getText()==null||this.singleQuantMenu_.getAmountOfIsotopes().getText().length()<1)
            canStartQuantification = false;
          else{
            try{
              Integer.parseInt(this.singleQuantMenu_.getAmountOfIsotopes().getText());
            }catch(NumberFormatException nfx){
              canStartQuantification = false;
            }
          }
          if (this.singleQuantMenu_.getAmountOfMatchingSearchIsotopes().getText()==null||this.singleQuantMenu_.getAmountOfMatchingSearchIsotopes().getText().length()<1)
            canStartQuantification = false;
          else{
            try{
              int isos = Integer.parseInt(this.singleQuantMenu_.getAmountOfMatchingSearchIsotopes().getText());
              if (isos<0)
                canStartQuantification = false;
            }catch(NumberFormatException nfx){
              canStartQuantification = false;
            }
          }
        }
        try{
          if (singleQuantMenu_.getCutoff().getText()!=null && singleQuantMenu_.getCutoff().getText().length()>0){
            cutoff = Float.parseFloat(singleQuantMenu_.getCutoff().getText().replaceAll(",", "."));
            LipidomicsConstants.getInstance().setRelativeMS1BasePeakCutoff(singleQuantMenu_.getCutoff().getText());
          }
        }catch(NumberFormatException ex){new WarningMessage(new JFrame(), "Error", "The cutoff value must be float format!"); canStartQuantification=false;}
        try{
          if (singleQuantMenu_.getRtShift().getText()!=null && singleQuantMenu_.getRtShift().getText().length()>0)
            rtShift = Float.parseFloat(singleQuantMenu_.getRtShift().getText().replaceAll(",", "."));
        }catch(NumberFormatException ex){new WarningMessage(new JFrame(), "Error", "The RT-shift value must be float format!"); canStartQuantification=false;}

//        if (this.singleQuantMenu_.getSearchUnknownTime().isSelected()){
//        }
      }      
      if (canStartQuantification){
        this.singleQuantMenu_.getQuantifyingPanel().setVisible(true);
        this.singleQuantMenu_.getStartQuantification().setEnabled(false);
        singleQuantMenu_.getSpinnerLabel().setVisible(true);
        String fileToTranslate = singleQuantMenu_.getSelectedMzxmlDirectory().getText();
        readFromRaw_ = false;
        boolean threadStarted = false;
        boolean aborted = false;
        int amountOfIsotopes = 0;
        int isotopesMustMatch = 0;
        if (this.singleQuantMenu_.getIsoValidation().isSelected()){
          amountOfIsotopes = Integer.parseInt(this.singleQuantMenu_.getAmountOfIsotopes().getText());
          isotopesMustMatch = Integer.parseInt(this.singleQuantMenu_.getAmountOfMatchingSearchIsotopes().getText());
        }
        
//        if (this.singleQuantMenu_.getSearchUnknownTime().isSelected()){
//          isotopesMustMatch = Integer.parseInt(this.singleQuantMenu_.getAmountOfMatchingSearchIsotopes().getText());
//        }
        String suffix = fileToTranslate.substring(fileToTranslate.lastIndexOf("."));
        if (suffix.equalsIgnoreCase(".RAW")||suffix.equalsIgnoreCase(".d")||suffix.equalsIgnoreCase(".wiff")){        
          File headerFile = new File(StringUtils.getChromFilePaths(fileToTranslate)[1]);
          if (!headerFile.exists()){
            File rawFile = new File(fileToTranslate);
            if ((rawFile.isFile()&& ((Settings.getReadWPath()!=null&&Settings.getReadWPath().length()>0)||(Settings.getMsConvertPath()!=null&&Settings.getMsConvertPath().length()>0)))||
              (rawFile.isDirectory() && ((suffix.equalsIgnoreCase(".RAW")&&((Settings.getMassWolfPath()!=null&&Settings.getMassWolfPath().length()>0)||((Settings.getMassPlusPlusPath()!=null&&Settings.getMassPlusPlusPath().length()>0))))
                                     ||   (suffix.equalsIgnoreCase(".d") &&Settings.getMsConvertPath()!=null&&Settings.getMsConvertPath().length()>0)))){
              String[] params = new String[3];
              boolean isMassPlusPlus = false;
              boolean watersMsConvert = false;
              if (rawFile.isFile()){
                if (!LipidomicsConstants.isMS2()||(Settings.getMsConvertPath()!=null&&Settings.getMsConvertPath().length()>0)){
                  if (Settings.getMsConvertPath()!=null&&Settings.getMsConvertPath().length()>0){
                    params = BatchQuantThread.getMsConvertParams(fileToTranslate);
                  }else if (Settings.getReadWPath()!=null&&Settings.getReadWPath().length()>0){
                    params[0] = Settings.getReadWPath();
                    params[1] = fileToTranslate;
                    params[2] = "p";
                  }
                }else{
                  aborted = true;
                  @SuppressWarnings("unused")
                  WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "For the MS/MS feature msconvert is required, because ReadW does not convert the MS/MS spectra of XCalibur correctly!");
                }  
              }
              if (rawFile.isDirectory()){
                if (suffix.equalsIgnoreCase(".RAW")){
                  if (!LipidomicsConstants.isMS2()||(Settings.getMassPlusPlusPath()!=null&&Settings.getMassPlusPlusPath().length()>0)||LipidomicsConstants.useMsconvertForWaters()){
                    String outputFile = fileToTranslate.substring(0,fileToTranslate.lastIndexOf("."))+"."+LipidomicsConstants.getIntermediateFileFormat();
                    if (LipidomicsConstants.useMsconvertForWaters()) {
                      params =BatchQuantThread.getMsConvertParamsWaters(fileToTranslate);
                      watersMsConvert = true;
                    }else if (Settings.getMassPlusPlusPath()!=null&&Settings.getMassPlusPlusPath().length()>0){
                      params = new String[8];
                      params[0] = Settings.getMassPlusPlusPath();
                      params[1] = "-in";
                      params[2] = fileToTranslate;
                      params[3] = "-out";
                      params[4] = LipidomicsConstants.getIntermediateFileFormat().toLowerCase(); //Not tested yet for mzML, also unsure whether it has to be lower case..
                      params[5] = outputFile;
                      params[6] = "-sample";
                      params[7] = "0";
                      if (LipidomicsConstants.isMS2()) isMassPlusPlus = true;
                    }else if (Settings.getMassWolfPath()!=null&&Settings.getMassWolfPath().length()>0){
                      params = new String[4];
                      params[0] = Settings.getMassWolfPath();
                      params[1] = "--"+LipidomicsConstants.getIntermediateFileFormat();
                      params[2] = fileToTranslate;
                      params[3] = outputFile;  
                    }
                  }else{
                    aborted = true;
                    @SuppressWarnings("unused")
                    WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "For the MS/MS feature Mass++ is required, because MassWolf does not convert the MS/MS spectra of MassLynx!");              
                  }
                }else if(suffix.equalsIgnoreCase(".d")){
                  if (Settings.getMsConvertPath()!=null&&Settings.getMsConvertPath().length()>0){
                    params =BatchQuantThread.getMsConvertParams(fileToTranslate);
                  }else{
                    aborted = true;
                    @SuppressWarnings("unused")
                    WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "For the conversion of MassHunter files \"msconvert\" is required!");
                  }                    
                }
              }
              if (params[0]!=null && params[0].length()>0){
                this.singleQuantMenu_.getProgressBar().setValue(5);
                this.singleQuantMenu_.getQuantifyingLabel().setText("Translating to "+LipidomicsConstants.getIntermediateFileFormat());
                rawmzThread_ = new RawToMzxmlThread(params,isMassPlusPlus,watersMsConvert);
                rawmzThread_.start();
                threadStarted = true;
              }  
            }
          }
        }
        if (!threadStarted && (fileToTranslate.endsWith(".mzXML") || fileToTranslate.endsWith(".mzML"))){
          File headerFile = new File(StringUtils.getChromFilePaths(fileToTranslate)[1]);
          if (!headerFile.exists()){
            this.singleQuantMenu_.getProgressBar().setValue(30);
            this.singleQuantMenu_.getQuantifyingLabel().setText("Translating to chrom");
            mzToChromThread_ = new MzxmlToChromThread(fileToTranslate,singleQuantMenu_.getNrProcessorsChrom());
            mzToChromThread_.start();
            threadStarted = true;
          }
        }
        if (!threadStarted && !aborted){
          boolean ionMode = false;
          if (this.singleQuantMenu_.getIonMode()!=null && ((String)singleQuantMenu_.getIonMode().getSelectedItem()).equalsIgnoreCase("+"))
            ionMode = true;
          this.singleQuantMenu_.getQuantifyingLabel().setText("Quantifying");
          this.singleQuantMenu_.getProgressBar().setValue(70);
          quantThread_ = new QuantificationThread(singleQuantMenu_.getSelectedMzxmlDirectory().getText(), singleQuantMenu_.getSelectedQuantDir().getText(),LipidDataAnalyzer.getResultFilePath(singleQuantMenu_.getSelectedMzxmlDirectory().getText(), singleQuantMenu_.getSelectedQuantDir().getText()),
              beforeTolerance,afterTolerance,amountOfIsotopes,isotopesMustMatch,
              this.singleQuantMenu_.getSearchUnknownTime().isSelected(),cutoff,rtShift,singleQuantMenu_.getNrProcessors(),ionMode,false);
          quantThread_.start();
          threadStarted = true;
        }else if (aborted){
          this.singleQuantMenu_.getQuantifyingPanel().setVisible(false);
          this.singleQuantMenu_.getStartQuantification().setEnabled(true);
          singleQuantMenu_.getSpinnerLabel().setVisible(false);          
        }
      }else{
//        boolean wasWarningMessage=false;
        if (singleQuantMenu_.getSelectedMzxmlDirectory().getText()==null||singleQuantMenu_.getSelectedMzxmlDirectory().getText().length()<1){
          @SuppressWarnings("unused")
          WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "You must specify a raw, mzXML, mzML or chrom file");
//          wasWarningMessage=true;
        }else if (singleQuantMenu_.getSelectedQuantDir().getText()==null||singleQuantMenu_.getSelectedQuantDir().getText().length()<1){
          @SuppressWarnings("unused")
          WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "You must specify a quant file");
//          wasWarningMessage=true;
        } else if (this.singleQuantMenu_.getIsoValidation().isSelected()){
          if (this.singleQuantMenu_.getAmountOfIsotopes().getText()==null||this.singleQuantMenu_.getAmountOfIsotopes().getText().length()<1){
            @SuppressWarnings("unused")
            WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "If you select the isotope option, please specify the amount of isotopes that you want to quantify (integer format)!");
//            wasWarningMessage=true;
          }else{
            try{
              Integer.parseInt(this.singleQuantMenu_.getAmountOfIsotopes().getText());
            }catch(NumberFormatException nfx){
              @SuppressWarnings("unused")
              WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "The number of isotopes must be integer format!");
//              wasWarningMessage=true;
            }
          }
          if (this.singleQuantMenu_.getAmountOfMatchingSearchIsotopes().getText()==null||this.singleQuantMenu_.getAmountOfMatchingSearchIsotopes().getText().length()<0){
            @SuppressWarnings("unused")
            WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "If you select the find molecules without retention-time option, please specify the amount of matching isotopes (integer format)!");
//            wasWarningMessage=true;
          }else{
            try{
              int isotopes = Integer.parseInt(this.singleQuantMenu_.getAmountOfMatchingSearchIsotopes().getText());
              if (isotopes<0){
                @SuppressWarnings("unused")
                WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "The number of matching isotopes must be >=0!");
//                wasWarningMessage=true;
              }
            }catch(NumberFormatException nfx){
              @SuppressWarnings("unused")
              WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "The number of matching isotopes must be integer format!");
//              wasWarningMessage=true;
            }
          }          
        }
//        if (!wasWarningMessage && this.singleQuantMenu_.getSearchUnknownTime().isSelected()){
//        }
      }
    }
    if (command.equalsIgnoreCase("startDisplay")){
    	if (selectedResultFile.getText()!=null)
        selectedResultFile.setText(selectedResultFile.getText().trim());
    	try
    	{
    		this.readResultFile(selectedResultFile.getText());
    	}
    	catch (ExcelInputFileException ex){}
    	startDisplay();
    }
    if (command.equalsIgnoreCase(DisplayTolerancePanel.UPDATE_QUANT_TOL_OF_CURRENTLY_SELECTED)){
      if (params_!=null){
        if (this.displaysMs2_){
          if (ms1ProbesWhileMs2Display_!=null){
            this.initANewViewer(params_,ms1ProbesWhileMs2Display_);
          }else
            this.initANewViewer(params_,getAllProbesFromParams(params_,null));
        }else{
          if (this.l2DPainter_!=null)
            this.initANewViewer(params_,l2DPainter_.getAllSelectedProbes());
          else
            this.initANewViewer(params_);
        }
      } else if (this.displayTolerancePanel_.getLockMzRange().isSelected()){
        this.currentSelectedSheet_ =  (String)selectedSheet_.getSelectedItem();
        initANewViewer(null);
      }
      lockRangeUpdateRequired_ = false;
    }
    if (command.equalsIgnoreCase("ZoomIn")||command.equalsIgnoreCase("ZoomMzIn")){
      float f1 = -1;
      float f2 = -1;
      Lipidomics2DPainter painter = null;
      if (command.equalsIgnoreCase("ZoomMzIn")){
        f1 = Float.parseFloat(mz_minTimeText_.getText());
        f2 = Float.parseFloat(mz_maxTimeText_.getText());
        painter = spectrumPainter_;
      }else{
        f1 = Float.parseFloat(m_minTimeText_.getText());
        f2 = Float.parseFloat(m_maxTimeText_.getText());
        painter = l2DPainter_;
      }
      if (f1 < 0 || f2 < 0 || f1 >= f2)
        return;
      if (painter.setMinDispTime2d(f1) == false)
        return;
      if (painter.setMaxDispTime2d(f2) == false)
        return;
      painter.repaint();
    }
    if (command.equalsIgnoreCase("SpectPlus")){
      selectCorrespondingSpectrum(true);
    }
    if (command.equalsIgnoreCase("SpectMinus")){
      selectCorrespondingSpectrum(false);
    }
    if (command.equalsIgnoreCase("SpectPlus")||command.equalsIgnoreCase("SpectMinus")){
      Float selectedRt = null;
      if (rtSelected_.getText()!=null && rtSelected_.getText().length()>0)
        selectedRt = Float.parseFloat(rtSelected_.getText());
      l2DPainter_.paintMs2Position(selectedRt);

    }
    if (command.equalsIgnoreCase("ZoomAll")||command.equalsIgnoreCase("ZoomMzAll")){
      Lipidomics2DPainter painter = null;
      if (command.equalsIgnoreCase("ZoomMzAll"))painter = spectrumPainter_;
      else painter = l2DPainter_;
      painter.zoomAll();
      painter.repaint();
    }
    if (command.equalsIgnoreCase("Dn")){
      l2DPainter_.previousChromatogram();
      l2DPainter_.repaint();
      viewer_.setCurrent2DMzRange(l2DPainter_.getLowerMz(), l2DPainter_.getUpperMz());
      viewer_.repaintColors();
    }
    if (command.equalsIgnoreCase("Up")){
      l2DPainter_.nextChromatogram();
      l2DPainter_.repaint();
      viewer_.setCurrent2DMzRange(l2DPainter_.getLowerMz(), l2DPainter_.getUpperMz());
      viewer_.repaintColors();
    }
    if (command.equalsIgnoreCase("Determine Area")) {
      l2DPainter_.determineArea(LipidomicsDefines.StandardValleyMethod);
      this.viewer_.setPaintableProbes(l2DPainter_.getAllSelectedProbes());
      this.viewer_.setTheShowSelectedWasOn(false);
      this.viewer_.repaintColors();
    }
    if (command.equalsIgnoreCase("Determine Area (Col)")) {
      l2DPainter_.determineArea(LipidomicsDefines.EnhancedValleyMethod);
      this.viewer_.setPaintableProbes(l2DPainter_.getAllSelectedProbes());
      this.viewer_.setTheShowSelectedWasOn(false);
      this.viewer_.repaintColors();
    }
    if (command.equalsIgnoreCase("Determine Area (Greedy)")) {
      l2DPainter_.determineArea(LipidomicsDefines.GreedySteepnessReductionMethod);
      this.viewer_.setPaintableProbes(l2DPainter_.getAllSelectedProbes());
      this.viewer_.setTheShowSelectedWasOn(false);
      this.viewer_.repaintColors();
    }
    if (command.equalsIgnoreCase("Determine Area (3D)")) {
      l2DPainter_.determineArea(LipidomicsDefines.Valley3DMethod);
      this.viewer_.setPaintableProbes(l2DPainter_.getAllSelectedProbes());
      this.viewer_.setTheShowSelectedWasOn(false);
      this.viewer_.repaintColors();
    }
    if (command.equalsIgnoreCase("Delete Area")) {
      l2DPainter_.deleteArea();
      this.viewer_.setPaintableProbes(l2DPainter_.getAllSelectedProbes());
      this.viewer_.setTheShowSelectedWasOn(false);
      this.viewer_.repaintColors();    
    }
    if (command.equalsIgnoreCase("acceptAreaSettings")){
      this.viewer_.setPaintableProbes(l2DPainter_.getAllSelectedProbes());
      this.viewer_.setTheShowSelectedWasOn(false);
      this.viewer_.repaintColors();      
    }
    if (command.equalsIgnoreCase("SaveSelectedAreas")){
      this.storeSelectedAreasToFile();
    }
    if (command.equalsIgnoreCase("addRemoveAbsoluteSettings")||subCommand.equalsIgnoreCase("removeSettings")){
      if (jButtonResultAbsQuant_.getText().equalsIgnoreCase("Add absolute settings") && command.equalsIgnoreCase("addRemoveAbsoluteSettings")){
        if (this.resultFiles_!=null&&this.resultFiles_.size()>0){
          ComparativeNameExtractor extractor;
          if (this.groupsPanel_.getGroups().size()>0){
            extractor = new ComparativeNameExtractor(this.resultFiles_,this.internalStandardSelection_.getText(),
                this.externalStandardSelection_.getText(),groupsPanel_.getGroups(),groupsPanel_.getGroupFiles());
          }else{
            extractor = new ComparativeNameExtractor(this.resultFiles_,this.internalStandardSelection_.getText(),
                this.externalStandardSelection_.getText());
          }
          
          try {
            extractor.parseInput(statisticsViewMode_, combineOxWithNonOx_);
            jButtonResultAbsQuant_.setText("Remove absolute settings");
            jButtonResultAbsQuant_.setToolTipText(TooltipTexts.STATISTICS_REMOVE_ABS_SETTINGS);
            resultWarningLabel_.setVisible(false);
            resultLoadButton_.setVisible(true); 
            resultSaveButton_.setVisible(true); 
            quantSettingsPanel_.showSettingsPanel(extractor);
          } catch (ExcelInputFileException e) {
            e.printStackTrace();
            new WarningMessage(new JFrame(), "Error", "Some of the input files contain invalid information!");
          } catch (LipidCombinameEncodingException e) {
            e.printStackTrace();
            new WarningMessage(new JFrame(), "Error", "There is something wrong with the lipid names in your result file!");
          }
        }else{
          @SuppressWarnings("unused")
          WarningMessage dlg = new WarningMessage(new JFrame(), "Error", "Please specify the files to analyze! There have to be at least one!!");
        }  

      }else{
        jButtonResultAbsQuant_.setText("Add absolute settings");
        jButtonResultAbsQuant_.setToolTipText(TooltipTexts.STATISTICS_ADD_ABS_SETTINGS);
        quantSettingsPanel_.hideSettingsPanel();
        resultWarningLabel_.setVisible(true);
        resultLoadButton_.setVisible(false); 
        resultSaveButton_.setVisible(false);
      }
    }
    if (command.equalsIgnoreCase("addRemoveCutoffSettings") || subCommand.equalsIgnoreCase("removeSettings")){
      if (jButtonResultCutoff_.getText().equalsIgnoreCase("Add cutoff settings") && command.equalsIgnoreCase("addRemoveCutoffSettings")){
        if (this.resultFiles_!=null&&this.resultFiles_.size()>0){
          ClassNamesExtractor extractor = new ClassNamesExtractor(this.resultFiles_);          
          try {
            extractor.parseInput(statisticsViewMode_, combineOxWithNonOx_);
            jButtonResultCutoff_.setText("Remove cutoff settings");
            jButtonResultCutoff_.setToolTipText(TooltipTexts.STATISTICS_REMOVE_CUTOFF_SETTINGS);
            cutoffWarningLabel_.setVisible(false);
            cutoffLoadButton_.setVisible(true); 
            cutoffSaveButton_.setVisible(true);
            cutoffSettingsPanel_.showSettingsPanel(extractor);
          }
          catch (ExcelInputFileException e) {
            e.printStackTrace();
            new WarningMessage(new JFrame(), "Error", "Some of the input files contain invalid information!");
          } catch (LipidCombinameEncodingException e) {
            e.printStackTrace();
            new WarningMessage(new JFrame(), "Error", "There is something wrong with the lipid names in your result file!");
          }
        }else{
          @SuppressWarnings("unused")
          WarningMessage dlg = new WarningMessage(new JFrame(), "Error", "Please specify the files to analyze! There have to be at least one!!");
        }  

      }else{
        jButtonResultCutoff_.setText("Add cutoff settings");
        jButtonResultCutoff_.setToolTipText(TooltipTexts.STATISTICS_ADD_CUTOFF_SETTINGS);
        cutoffSettingsPanel_.hideSettingsPanel();
        cutoffWarningLabel_.setVisible(true);
        cutoffLoadButton_.setVisible(false); 
        cutoffSaveButton_.setVisible(false);
      }
    }
    if (command.equalsIgnoreCase("saveCutoffSettings")){
      int returnVal = saveCutoffSettingsFileChooser_.showSaveDialog(new JFrame());
      if (returnVal == JFileChooser.APPROVE_OPTION) {
        File fileToStore = saveCutoffSettingsFileChooser_.getSelectedFile();
        boolean store = true;
        if (fileToStore.exists()){
          if (JOptionPane.showConfirmDialog(this, "The file "+fileToStore.getName()+" exists! Replace existing file?") != JOptionPane.YES_OPTION)
            store = false;
        }else{
          if (fileToStore.getName().indexOf(".")==-1){
            fileToStore = new File(fileToStore.getAbsoluteFile()+".cut.xml");
            if (fileToStore.exists()){
              if (JOptionPane.showConfirmDialog(this, "The file "+fileToStore.getName()+" exists! Replace existing file?") != JOptionPane.YES_OPTION)
                store = false;
            }  
          }  
        }    
        if (store){
          Exception ex = null;
          try {
            CutoffSettingsWriter writer = new CutoffSettingsWriter(fileToStore.getAbsolutePath(), this.cutoffSettingsPanel_.getMaxIsotope(),  cutoffSettingsPanel_.getCutoffsAsString());
            writer.writeXMLFile();
          }
          catch (ParserConfigurationException e) {
            ex = e;
          } catch (TransformerConfigurationException e) {
            ex = e;
          }
          catch (TransformerException e) {
            ex = e;
          }
          catch (IOException e) {
            ex = e;
          } catch (AbsoluteSettingsInputException e){
            ex = e;            
          }
          if (ex!=null){
            ex.printStackTrace();
            new WarningMessage(new JFrame(), "Error", ex.getMessage());
          }

        }
      } 
      return;
    }
    if (command.equalsIgnoreCase("saveAbsoluteSettings")){
      int returnVal = saveAbsSettingsFileChooser_.showSaveDialog(new JFrame());
      if (returnVal == JFileChooser.APPROVE_OPTION) {
        File fileToStore = saveAbsSettingsFileChooser_.getSelectedFile();
        boolean store = true;
        if (fileToStore.exists()){
          if (JOptionPane.showConfirmDialog(this, "The file "+fileToStore.getName()+" exists! Replace existing file?") != JOptionPane.YES_OPTION)
            store = false;
        }else{
          if (fileToStore.getName().indexOf(".")==-1){
            fileToStore = new File(fileToStore.getAbsoluteFile()+".wqs.xml");
            if (fileToStore.exists()){
              if (JOptionPane.showConfirmDialog(this, "The file "+fileToStore.getName()+" exists! Replace existing file?") != JOptionPane.YES_OPTION)
                store = false;
            }  
          }  
        }    
        if (store){
          AbsoluteQuantSettingsWholeWriter writer = new AbsoluteQuantSettingsWholeWriter(fileToStore.getAbsolutePath(),this.quantSettingsPanel_);
          Exception ex = null;
          try {
              writer.writeXMLFile();
          }
          catch (ParserConfigurationException e) {
            ex = e;
          } catch (TransformerConfigurationException e) {
            ex = e;
          }
          catch (TransformerException e) {
            ex = e;
          }
          catch (IOException e) {
            ex = e;
          }
          if (ex!=null){
            ex.printStackTrace();
            new WarningMessage(new JFrame(), "Error", ex.getMessage());
          }

        }
         //this.selectedResultFile.setText(text);

      } 
      return;
    }
    if (command.equalsIgnoreCase("loadCutoffSettings")){
      new WarningMessage(new JFrame(), "Warning", "Pay attention that the loading of the file will erase your previous \"cutoff-settings\"");
      if (loadCutoffSettingsFileChooser_.showOpenDialog(this) == JFileChooser.APPROVE_OPTION){
        File file = loadCutoffSettingsFileChooser_.getSelectedFile();
        if (file.exists()){
          try {
            CutoffSettingsReader reader = new CutoffSettingsReader(file.getAbsolutePath(),this.cutoffSettingsPanel_);
            reader.parseXMLFile();
          }
          catch (Exception e) {
            e.printStackTrace();
            new WarningMessage(new JFrame(), "Error", e.getMessage());
          }
        }
      }
      return;
    }

    if (command.equalsIgnoreCase("loadAbsoluteSettings")){
      new WarningMessage(new JFrame(), "Warning", "Pay attention that the loading of the file will erase your previous \"absolute quant-settings\"");
      if (loadAbsSettingsFileChooser_.showOpenDialog(this) == JFileChooser.APPROVE_OPTION){
        File file = loadAbsSettingsFileChooser_.getSelectedFile();
        if (file.exists()){
          try {
            AbsoluteQuantSettingsWholeReader reader = new AbsoluteQuantSettingsWholeReader(file.getAbsolutePath(),this.quantSettingsPanel_);
            reader.parseXMLFile();
          }
          catch (Exception e) {
            e.printStackTrace();
            new WarningMessage(new JFrame(), "Error", e.getMessage());
          }
        }
      }
      return;
    }
    if (command.equalsIgnoreCase("ApplyOtherMachineSettings")){
      try {
        Settings.applySettings((String)msMachineTypes_.getSelectedItem());
        Settings.applyFragmentationSettings((String)fragmentationSettings1_.getSelectedItem(),(String)fragmentationSettings2_.getSelectedItem());
        frame_.setTitle(getFrameTitleString());
      }
      catch (SettingsException e) {
        e.printStackTrace();
        new WarningMessage(new JFrame(), "Error", e.getMessage());
      }      
      return;
    }
    if (command.equalsIgnoreCase("SaveOtherMachineSettings")){
      try {
        Settings.saveMachineSettings((String)msMachineTypes_.getSelectedItem());
        Settings.saveFragmentationSettings((String)fragmentationSettings1_.getSelectedItem(),(String)fragmentationSettings2_.getSelectedItem());
        frame_.setTitle(getFrameTitleString());
      }
      catch (SettingsException e) {
        e.printStackTrace();
        new WarningMessage(new JFrame(), "Error", e.getMessage());
      }      
      return;
    }
    if (command.equalsIgnoreCase(CHANGE_SEPARATE_RT_STATUS)){
    	if (LipidomicsConstants.isShotgun()==1&&separateHitsByRT_.isSelected())
    	{
    		new WarningMessage(new JFrame(), "Warning", "Shotgun settings are selected. RT grouping is not enabled for shotgun data!");
    		separateHitsByRT_.setSelected(false);
    	}
    	else
    	{
    		rtGroupingTime_.setEnabled(separateHitsByRT_.isSelected());
        rtTimeUnit_.setEnabled(separateHitsByRT_.isSelected());
    	}
    } else if  (command.equalsIgnoreCase(ExportPanel.EXPORT_PNG)){
      exportFileChooser_.setFileFilter(new FileNameExtensionFilter("PNG (*.png)","png"));
      int returnVal = exportFileChooser_.showSaveDialog(new JFrame());
      if (returnVal == JFileChooser.APPROVE_OPTION) {
        File fileToStore = exportFileChooser_.getSelectedFile();
        @SuppressWarnings("rawtypes")
        Vector results = StaticUtils.checkFileStorage(fileToStore,"png",this);
        fileToStore = (File)results.get(0);
        if ((Boolean)results.get(1)){
          try {
            ImageIO.write(spectrumPainter_.getImage(), "PNG", fileToStore);
          }catch (IOException ex){ new WarningMessage(new JFrame(), "Error", ex.getMessage());}
        }
      }
    }else if (command.equalsIgnoreCase(ExportPanel.EXPORT_SVG)){
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
            if (!exportChromatogramsFromDRView_)
            {
            	spectrumPainter_.draw2DDiagram(svgGenerator);
            }
            else
            {
            	l2DPainter_.draw2DDiagram(svgGenerator);
            }
            boolean useCSS = true; // we want to use CSS style attributes
            BufferedOutputStream stream = new BufferedOutputStream(new FileOutputStream(fileToStore));
            Writer out = new OutputStreamWriter(stream, "UTF-8");
            svgGenerator.stream(out, useCSS);
            out.close();
          }catch (IOException ex){ new WarningMessage(new JFrame(), "Error", ex.getMessage());}
        }  
      }      
    }else if (command.equalsIgnoreCase("mgf")){
      exportFileChooser_.setFileFilter(new FileNameExtensionFilter("MGF (*.mgf)","mgf"));
      int returnVal = exportFileChooser_.showSaveDialog(new JFrame());
      if (returnVal == JFileChooser.APPROVE_OPTION) {
        File fileToStore = exportFileChooser_.getSelectedFile();
        @SuppressWarnings("rawtypes")
        Vector results = StaticUtils.checkFileStorage(fileToStore,"mgf",this);
        fileToStore = (File)results.get(0);
        if ((Boolean)results.get(1)){
          try {
            String chromFile = StringUtils.getJustFileName(selectedChromFile.getText());
            chromFile = chromFile.substring(0,chromFile.length()-".chrom".length());
            ((Lipidomics2DSpectraChromPainter)spectrumPainter_).writeMgf(fileToStore,StringUtils.getJustFileName(chromFile));
          }catch (IOException ex){ new WarningMessage(new JFrame(), "Error", ex.getMessage());}
        }  
      }      
    }else if (command.equalsIgnoreCase("AcceptMSnRecalculation")){
      updateLipidParameterSet(recalcDialog_.getResult(), currentSelected_);
      recalcDialog_ = null;
    }else if (command.equalsIgnoreCase("DeclineMSnRecalculation")){
      recalcDialog_ = null;
    }else if (command.equalsIgnoreCase("AcceptChangedRT")){
      LipidParameterSet set = getAnalyteInTableAtPosition(currentSelected_);
      LipidParameterSet originalSet = new LipidParameterSet(set);
      set.setPreciseRt(editRtDialog_.getRt());
      try {
        storeResultsToExcel(currentSelected_);
      } catch (ExportException ex) {
        set = originalSet;
      }
      editRtDialog_ = null;
    }else if (command.equalsIgnoreCase("DeclineChangedRT")){
      editRtDialog_ = null;
    } else if (command.equals(EditOmegaAssignmentJTable.getSaveChangesCommand())) {
      storeOmegaResultsToExcel(currentSelected_);
    } else if (command.equalsIgnoreCase("AcceptExportSettings")){
      exportSettings_.setVisible(false);
      if (exportSettingsGroup_!=null)
        exportSettingsGroup_.setVisible(false);
    }
    if (command.equalsIgnoreCase("Start AI"))
    {
//    	String suffix = "_shortened.xlsx";
    	String suffix = ".xlsx";
    	File file = new File("D:\\Collaborator_Files\\Kathi\\Paper3\\LDA_extension\\Description\\Gangliosides_targets\\Target_list_gangliosides_adducts"+suffix);
    	String outPath = "D:\\Collaborator_Files\\Kathi\\Paper3\\LDA_extension\\Description\\Gangliosides_targets\\TestChrom\\spectra"+suffix;
    	String outPathMerged = "D:\\Collaborator_Files\\Kathi\\Paper3\\LDA_extension\\Description\\Gangliosides_targets\\TestChrom\\spectra_merged"+suffix;
    	ExcelTargetListParser parser = new ExcelTargetListParser(file);
    	SpectraIdentifier identifier = new SpectraIdentifier("D:\\Collaborator_Files\\Kathi\\Paper3\\LDA_extension\\data\\LCMS_STD_Glycolipids_data");
    	try
    	{
    		parser.parse();
    		ArrayList<SpectrumContainer> spectra = identifier.identifySpectra(parser);
    		SpectraTextExporter exporter = new SpectraTextExporter(spectra, outPath);
    		exporter.exportSpectra();
    		SpectraInterpreter interpreter = new SpectraInterpreter(spectra);
    		ArrayList<CombinedSpectrumContainer> mergedSpectra = interpreter.interpretSpectra();
    		exporter = new SpectraTextExporter(outPathMerged, mergedSpectra);
    		exporter.exportSpectra();
    	}
    	catch (IOException | QuantificationException ex)
    	{
    		ex.printStackTrace();
    	}
    }
  }
  
  private void startDisplay()
  {
  	if (selectedChromFile.getText()!=null)
      selectedChromFile.setText(selectedChromFile.getText().trim());
    if (selectedChromFile.getText()!=null&&selectedChromFile.getText().length()>0 &&
        selectedResultFile.getText()!=null&&selectedResultFile.getText().length()>0){
      String[] chromPaths = StringUtils.getChromFilePaths(selectedChromFile.getText());
      String pureFile = chromPaths[0].substring(0,selectedChromFile.getText().lastIndexOf("."));
      if (StaticUtils.existChromFiles(pureFile) && StaticUtils.existsFile(selectedResultFile.getText())){
        try {
          reader_ = new ChromatogramReader(chromPaths[1], chromPaths[2], chromPaths[3],  chromPaths[0],LipidomicsConstants.isSparseData(),LipidomicsConstants.getChromSmoothRange());
          analyzer_ = new LipidomicsAnalyzer(chromPaths[1], chromPaths[2], chromPaths[3],  chromPaths[0],false);           
          
          currentSelected_ = -1;
          currentSelectedSheet_ = "";
          QuantificationThread.setAnalyzerProperties(analyzer_);
          analyzer_.setGeneralBasePeakCutoff(0f);
          
          Component[] components = tablePanel_.getComponents();
          for (Component component : components) {
            if (component.equals(displayTolerancePanel_)) {
              tablePanel_.remove(component);
            }
          }
          boolean showOmegaEditingTools = LipidParameterSet.isOmegaInformationAvailable() || Settings.getAlwaysEditOmega();
          displayTolerancePanel_ = new DisplayTolerancePanel(this,showOmegaEditingTools);
          tablePanel_.add(displayTolerancePanel_,BorderLayout.SOUTH);
          this.updateSheetSelectionList();
          this.params_ = null;
          RuleDefinitionInterface.clearCacheDir();
          lockRangeUpdateRequired_ = true;
          displayTolerancePanel_.getShowMSnEvidence().setSelected(false);
          displayTolerancePanel_.getShowChainEvidence().setSelected(false);
          displayTolerancePanel_.getShowChainEvidence().setEnabled(false);
        }
        catch (CgException e) {
          @SuppressWarnings("unused")
          WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "The chrom file you specified is not valid");
          e.printStackTrace();
        }
      } else {
        if (!StaticUtils.existsFile(selectedResultFile.getText()))
          new WarningMessage(new JFrame(), "Warning", "The results file \""+StaticUtils.extractFileName(selectedResultFile.getText())+"\" does not exist!");
      }
    }else{
      if (selectedChromFile.getText()==null||selectedChromFile.getText().length()<1){
        @SuppressWarnings("unused")
        WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "You must specify a chrom file");
      }else if (selectedResultFile.getText()==null||selectedResultFile.getText().length()<1){
        @SuppressWarnings("unused")
        WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "You must specify a result file");
      }
    }
  }
  
  public void showExportSettingsDialog(boolean grouped){
    if (grouped)
      exportSettingsGroup_.setVisible(true);
    else
      exportSettings_.setVisible(true);
  }
  
  /**
   * Adjusts display settings of all heatmaps to a reference object as far as applicable.
   */
  public void applySettingsToAllClasses(ResultDisplaySettings settings)
  {
    for (String molGroup : analysisModule_.getResults().keySet())
    {
    	if (heatmaps_ != null)
    	{
    		HeatMapDrawing drawing = heatmaps_.get(molGroup);
    		drawing.adjustDisplaySettings(settings);
    	}
    	
    	if (groupHeatmaps_ != null)
    	{
    		HeatMapDrawing drawing = groupHeatmaps_.get(molGroup);
    		if (drawing != null)
    		{
    			drawing.adjustDisplaySettings(settings);
    		}
    	}
    }
  }
  
  /**
   * updates and saves a LipidParameterSet object at a certain position in the displayed table
   * @param newOne the new LipidParameterSet object
   * @param position the position in the table where the object shall be updated
   */
  private void updateLipidParameterSet(LipidParameterSet newOne, int position){
    QuantificationResult originalResult = new QuantificationResult(result_);
    int originalPosition = resultPositionToOriginalLoopkup_.get(currentSelected_);
    result_.getIdentifications().get(currentSelectedSheet_).remove(originalPosition);
    result_.getIdentifications().get(currentSelectedSheet_).add(originalPosition,newOne);
    try {
      storeResultsToExcel(currentSelected_);
    } catch (ExportException ex) {
      result_ = originalResult;
    }
  }
  
  /**
   * stores the current results back to the Excel file, only if successful the selection table is updated
   * @param position position to change the selection to
   */
  private void storeResultsToExcel(int position) throws ExportException {
    try {
      QuantificationResultExporter.writeResultsToExcel(selectedResultFile.getText(), result_);
      this.updateResultListSelectionTable();
      this.displayTable_.changeSelection(position, 1, false, false);
      listSelectionChanged(position);
    } catch (ExportException ex) {
      new WarningMessage(new JFrame(), "Error", ex.getMessage());
      ex.printStackTrace();
      throw new ExportException(ex);
    }   
  }
  
  private void storeOmegaResultsToExcel(int position) {
    try {
      QuantificationResultExporter.writeResultsToExcel(selectedResultFile.getText(), result_);
      EditOmegaAssignmentJTable.setSaveLipidParameterSet(true);
      updateResultListSelectionTable();
      displayTable_.changeSelection(position, 1, false, false);
      listSelectionChanged(position);
    } catch (ExportException ex) {
      EditOmegaAssignmentJTable.setSaveLipidParameterSet(false);
      new WarningMessage(new JFrame(), "Error", ex.getMessage());
      ex.printStackTrace();
    }
  }
  
  private void selectCorrespondingSpectrum(boolean next){
    RuleDefinitionInterface rdi = null;
    if (topSplitPane_.getTopComponent() instanceof RuleDefinitionInterface) rdi = (RuleDefinitionInterface) topSplitPane_.getTopComponent();
    LipidParameterSet param = null;
    Hashtable<Integer,Vector<RangeColor>> rangeColors = new Hashtable<Integer,Vector<RangeColor>>();
    if (rdi != null) {
      int specNumber = getSelectedSpectrumNumber();
      if (next) specNumber++;
      else specNumber--;
      try {
       param = rdi.testForMSnDetection(spectrumPainter_.getMs2LevelSpectrumSelected(specNumber));
       rangeColors = StaticUtils.createRangeColorVOs(param, ((LipidomicsTableModel)displayTable_.getModel()).getMSnIdentificationName(currentSelected_),
           result_.getFaHydroxyEncoding(), result_.getLcbHydroxyEncoding(), areTheseAlex123MsnFragments());
      }
      catch (NoRuleException | IOException
          | SpectrummillParserException | CgException | ChemicalFormulaException | LipidCombinameEncodingException e) {
        e.printStackTrace();
      }
      //if it is a rules exception, then a new rule is currently entered - no need to report an error
      catch (RulesException  | HydroxylationEncodingException e) {
      }
   }
    if (next){
      if (param!=null) spectrumPainter_.nextSpectrum(param,rangeColors);
      else spectrumPainter_.nextSpectrum();
    } else {
      if (param!=null) spectrumPainter_.previousSpectrum(param,rangeColors);
      else spectrumPainter_.previousSpectrum();
    }
    this.spectrumSelected_.setText(spectrumPainter_.getSpectSelectedText());
    this.rtSelected_.setText(spectrumPainter_.getRtSelectedText());
    displayPrecursorMasses(spectrumPainter_.getPrecursorMassSelected());
    this.msLevelSelected_.setText(spectrumPainter_.getMsLevelSelected());
    if (rdi==null){
      float[] range = spectrumPainter_.getRTRange();
      viewer_.setCurrent2DTimeRange(range[0], range[1]);
      try{viewer_.repaintColors();}catch(Exception ex){}
    }
  }
  
  @SuppressWarnings("unchecked")
  private void acceptResultFiles() throws ExcelInputFileException
  {
//    long timeMillis_0 = System.currentTimeMillis();
//    System.out.println("fastExcel Start!");
    expDisplayNamesLookup_ = new Hashtable<String,String>();
    if (this.resultFiles_!=null&&this.resultFiles_.size()>0){
      AbsoluteSettingsVO absSettingVO = null;
      Hashtable<String,Double> cutoffValues = null;
      int maxCutoffIsotope = -1;
      if (jButtonResultAbsQuant_.getText().equalsIgnoreCase("Remove absolute settings")){
        try {
          absSettingVO = quantSettingsPanel_.getSettingsVO();
        }
        catch (AbsoluteSettingsInputException e) {
          new WarningMessage(new JFrame(), "Error", e.getMessage());
          return;
        }
      }
      if (jButtonResultCutoff_.getText().equalsIgnoreCase("Remove cutoff settings")){
        try {
          cutoffValues = cutoffSettingsPanel_.getCutoffs();
          maxCutoffIsotope = cutoffSettingsPanel_.getMaxIsotope();
        }
        catch (AbsoluteSettingsInputException e) {
          new WarningMessage(new JFrame(), "Error", e.getMessage());
          return;
        }
        
      }
      LinkedHashMap<String,Integer> classSequence = null;
      LinkedHashMap<String,Vector<String>> correctAnalyteSequence = null;
      Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects = null;
      if (correctOrderFile_.getText()!=null && correctOrderFile_.getText().length()>0){
        try {
          boolean ionMode = false;
          if (this.ionModeOrder_!=null && ((String)ionModeOrder_.getSelectedItem()).equalsIgnoreCase("+"))
            ionMode = true;
          @SuppressWarnings("rawtypes")
          Vector quantInfo = QuantificationThread.getCorrectAnalyteSequence(correctOrderFile_.getText(),ionMode);
          classSequence = (LinkedHashMap<String,Integer>)quantInfo.get(0);
          correctAnalyteSequence = (LinkedHashMap<String,Vector<String>>)quantInfo.get(1);
          quantObjects = (Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>>)quantInfo.get(4);
          
        }
        catch (Exception e) {
          new WarningMessage(new JFrame(), "Error", e.getMessage());
        }
      }
      separateHitsByRT_.setToolTipText(TooltipTexts.STATISTICS_SEPARATE_RT);
      double expRtGroupingTime = -1d;
      
      if (separateHitsByRT_.isSelected()) expRtGroupingTime = Double.valueOf(rtGroupingTime_.getText());

      if (this.groupsPanel_.getGroups().size()>0){
        analysisModule_ = new ComparativeAnalysis(this.resultFiles_,this.internalStandardSelection_.getText(),
            this.externalStandardSelection_.getText(), absSettingVO, cutoffValues, maxCutoffIsotope, groupsPanel_.getGroups(),
            groupsPanel_.getGroupFiles(),classSequence,correctAnalyteSequence,quantObjects,expRtGroupingTime);
      }else{
        analysisModule_ = new ComparativeAnalysis(this.resultFiles_,this.internalStandardSelection_.getText(),
            this.externalStandardSelection_.getText(), absSettingVO, cutoffValues, maxCutoffIsotope, 
            classSequence,correctAnalyteSequence,quantObjects,expRtGroupingTime);
      }
      try {
//        long before = System.currentTimeMillis();
        analysisModule_.parseInput(statisticsViewMode_,combineOxWithNonOx_);
//        System.out.println(String.format("Total time required by the reading: %s !", 
//            (System.currentTimeMillis()-before)/1000.0));
        analysisModule_.calculateStatistics();
        analysisModule_.getNrOfChainsOfClass();
        this.expDisplayNamesLookup_ = analysisModule_.getExpNames();
        if (this.groupsPanel_.getGroups().size()>0){
          groupDisplayNamesLookup_ = new Hashtable<String,String>();
          for (String groupName : this.groupsPanel_.getGroups()){
            groupDisplayNamesLookup_.put(groupName, groupName);
          }
        }
      }
      catch (ExcelInputFileException e) {
      	throw new ExcelInputFileException("Some of the input files contain invalid information!");
      }
      catch (Exception e) {
      	throw new ExcelInputFileException("An Excel file returns the following failure: "+e.getMessage());
      }
      this.generateHeatMaps();
    }else{
      new WarningMessage(new JFrame(), "Error", "Please specify the files to analyze!");
    }
//    long timeMillis_1 = System.currentTimeMillis();
//    System.out.println(String.format("Total time required by fastExcel: %s !", 
//        (timeMillis_1-timeMillis_0)/1000.0));
  }
  
  private void removeResultTabComponentsExceptFirst()
  {
  	while (resultTabs_.getTabCount() > 1) 
  	{
  		resultTabs_.remove(1);
  	}
  }
  
  private void generateHeatMaps(){
    heatmaps_ = new Hashtable<String,HeatMapDrawing>();
    Hashtable<String,Hashtable<String,Hashtable<String,ResultCompVO>>> analysisResults = analysisModule_.getResults();
    Hashtable<String,Hashtable<String,Integer>> corrTypeISLookup = analysisModule_.getCorrectionTypeISLookup();
    Hashtable<String,Hashtable<String,Integer>> corrTypeESLookup = analysisModule_.getCorrectionTypeESLookup();
    Vector<String> expNames = analysisModule_.getExpNamesInSequence();
    molBarCharts_ = new Hashtable<String,JTabbedPane>();
    
    // for the groupedValues
    Hashtable<String,Hashtable<String,Hashtable<String,ResultCompVO>>> groupedAnalysisResults = analysisModule_.getGroupedResults();
    groupHeatmaps_ = new Hashtable<String,HeatMapDrawing>();
    Hashtable<String,ResultDisplaySettings> displaySettingHash = new Hashtable<String,ResultDisplaySettings>();
    
    colorChooserDialog_ = new ColorChooserDialog(new JFrame(),"",expNames,groupsPanel_.getGroups(),this);
    exportSettings_ = new ExportSettingsPanel(false,this);
    exportSettingsGroup_ = null;
    if (this.groupsPanel_.getGroups().size()>0)
      exportSettingsGroup_ = new ExportSettingsPanel(true,this);
    
//    long before = System.currentTimeMillis();
    
    ExecutorService threadpool = Executors.newFixedThreadPool(Math.min(analysisResults.keySet().size(), getAmountOfProcessorsPreferred()));
    //JTabbedPane is not threadsafe. Thus, the panels are stored in this collection to be added iteratively.
    Hashtable<String,JPanel> jPanels = new Hashtable<String,JPanel>();
    for (String molGroup : analysisResults.keySet())
    {
//    	HeatMapBuilder builder = new HeatMapBuilder(displaySettingHash, molGroup, jPanels);
//    	builder.run();
    	threadpool.execute(new HeatMapBuilder(displaySettingHash, molGroup, jPanels));
    }
    threadpool.shutdown();
    try { threadpool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS); } catch (InterruptedException e) {}
    
    for (String molGroup : heatmaps_.keySet())
    {
    	resultTabs_.addTab(molGroup,jPanels.get(molGroup));
      resultTabs_.setToolTipTextAt(resultTabs_.indexOfTab(molGroup), TooltipTexts.TABS_RESULTS_GROUP+molGroup+"</html>");
    }
    
//    System.out.println(String.format("Time required by the heatmap gen in total: %s !", 
//    		(System.currentTimeMillis()-before)/1000.0));
    
    if (analysisResults.size()>0)
      resultTabs_.setSelectedIndex(1);
    if (analysisResults.size()>1){
    	String overview = "Overview";
      classOverviewPanel_ = new ClassesOverviewPanel(expNames,this,analysisModule_,displaySettingHash,heatmaps_,corrTypeISLookup,corrTypeESLookup,analysisResults,colorChooserDialog_);
      resultTabs_.addTab(overview, classOverviewPanel_);
      resultTabs_.setToolTipTextAt(resultTabs_.indexOfTab(overview), TooltipTexts.TABS_RESULTS_OVERVIEW);
      if (this.groupsPanel_.getGroups().size()>0){
      	String overviewGroup = "O.view-Group";
        classOverviewGroupPanel_ = new ClassesOverviewPanel(groupsPanel_.getGroups(),this,analysisModule_,displaySettingHash,heatmaps_,corrTypeISLookup,corrTypeESLookup,groupedAnalysisResults,colorChooserDialog_);
        resultTabs_.addTab(overviewGroup, classOverviewGroupPanel_);
        resultTabs_.setToolTipTextAt(resultTabs_.indexOfTab(overviewGroup), TooltipTexts.TABS_RESULTS_OVERVIEW_GROUPS);
      }
    }

  }
  
  private void readResultFile(String filePath) throws ExcelInputFileException{
    readResultFile(filePath,false);
  }
  
  private void readResultFile(String filePath,boolean keepOrder) throws ExcelInputFileException{
    resultsShowModification_ = new Hashtable<String,Boolean>();
    if (!keepOrder) orderResultsType_ = new Hashtable<String,Integer>();
    result_ = LDAResultReader.readResultFile(filePath,  resultsShowModification_);
    
    //start: added via the oxidized lipids extension
    originalResult_ = new QuantificationResult(result_); //changed to deep copy instead of reading the file again => TODO: make sure this is fine
    
    mSnResult_ = LDAResultReader.readResultFile(filePath,  resultsShowModification_);
    chainResult_ = LDAResultReader.readResultFile(filePath,  resultsShowModification_);
    Hashtable<String,Vector<LipidParameterSet>> mSnHash = new Hashtable<String,Vector<LipidParameterSet>>();
    Hashtable<String,Vector<LipidParameterSet>> chainHash = new Hashtable<String,Vector<LipidParameterSet>>();
    
    for (String lipidClass : result_.getIdentifications().keySet()) {
    	Vector<LipidParameterSet> params = result_.getIdentifications().get(lipidClass);
    	Vector<LipidParameterSet> mSnSets = new Vector<LipidParameterSet>();
    	Vector<LipidParameterSet> chainSets = new Vector<LipidParameterSet>();
    	for (LipidParameterSet param  : params){
  			if (param instanceof LipidomicsMSnSet) {
  				mSnSets.add(param);
  				LipidomicsMSnSet msn_param = (LipidomicsMSnSet) param;
  				if(!msn_param.getChainFragments().isEmpty())
  				{
  					chainSets.add(param);
  				}
  			}
  		}
    	mSnHash.put(lipidClass, mSnSets);
    	chainHash.put(lipidClass, chainSets);
    }
    mSnResult_.setIdentifications(mSnHash);
    chainResult_.setIdentifications(chainHash);
    //end: added via the oxidized lipids extension
    
    if (result_.getConstants()!=null && result_.getConstants().getShotgun()>LipidomicsConstants.SHOTGUN_FALSE)
      disableChromatographyFeatures();
    else
      enableChromatographyFeatures();
  }
  
  public void initANewViewer(int position){
    this.currentSelected_ = position;
    this.initANewViewer(getAnalyteInTableAtPosition(currentSelected_));
  }

  
  private void initANewViewer(LipidParameterSet params){
    this.initANewViewer(params,null);
  }
  
  private void initANewViewer(LipidParameterSet params,Vector<Vector<CgProbe>> previouslyselectedProbes){
    try {
      displaysMs2_ = false;
      int charge = 1;
      float startFloat;
      float stopFloat;
      float currentIsotopicMass = 0;
      if (params!=null) {
        if (params.ProbeCount()>0)
          charge = params.Probe(0).Charge;
        else if (params.getCharge()!=null&&params.getCharge()>1)
          charge = params.getCharge();
        currentIsotopicMass = params.Mz[0]+(LipidomicsConstants.getNeutronMass()*Integer.parseInt((String)this.isotope_.getSelectedItem())/(float)charge);
      }
      
      if (this.displayTolerancePanel_.getLockMzRange().isSelected()) {
        startFloat = Float.parseFloat(displayTolerancePanel_.getDisplayMzStart().getText());
        stopFloat = Float.parseFloat(displayTolerancePanel_.getDisplayMzStop().getText());
        if (startFloat>stopFloat) {
          new WarningMessage(new JFrame(), "Error", "The stop value of the \""+DisplayTolerancePanel.DISPLAY_LOCK_MZ_TEXT+"\" cannot be smaller than the start m/z value");
          return;
        }
      } else {
        startFloat = currentIsotopicMass-Float.parseFloat(this.displayTolerancePanel_.getDisplayMinusTolerance().getText());
        stopFloat = currentIsotopicMass+Float.parseFloat(this.displayTolerancePanel_.getDisplayPlusTolerance().getText());
      }
      float startRt = 0f;
      if (this.displayTolerancePanel_.getDisplayRtStart().getText()!=null && this.displayTolerancePanel_.getDisplayRtStart().getText().length()>0)
        startRt = Float.parseFloat(this.displayTolerancePanel_.getDisplayRtStart().getText());
      startRt = 60f*startRt;
      float stopRt = 0f;
      if (this.displayTolerancePanel_.getDisplayRtStop().getText()!=null && this.displayTolerancePanel_.getDisplayRtStop().getText().length()>0)
        stopRt = Float.parseFloat(this.displayTolerancePanel_.getDisplayRtStop().getText());
      stopRt = 60f*stopRt;
      
      String[] rawLines = reader_.getRawLines(startFloat, stopFloat, result_.getMsLevels().get(currentSelectedSheet_));
      Hashtable<Integer,Float> rtTimes = reader_.getRetentionTimesOriginal();

      MSMapViewer viewer = MSMapViewerFactory.getMSMapViewer(rawLines, rtTimes,startFloat,stopFloat,
          startRt,stopRt,reader_.getMultiplicationFactorForInt_()/reader_.getLowestResolution_(),5f,this,
          MSMapViewer.DISPLAY_TIME_MINUTES,false);
      viewer.setViewerSettings(true, true, true,LipidomicsConstants.getThreeDViewerDefaultMZResolution(),LipidomicsConstants.getThreeDViewerDefaultTimeResolution());
      //writeDisplayDataToExcelFormat(rawLines, rtTimes, startFloat,stopFloat);
      
      Vector<CgProbe> storedProbes = new Vector<CgProbe>();
      Vector<CgProbe> selectedProbes = new Vector<CgProbe>();
      if (params!=null) {
        Vector<Vector<CgProbe>> allProbes = getAllProbesFromParams(params,previouslyselectedProbes);
        storedProbes = allProbes.get(0);
        selectedProbes = allProbes.get(1);
        viewer.setPaintableProbes(allProbes);
        viewer.setCurrent2DMzRange(currentIsotopicMass-params.LowerMzBand, currentIsotopicMass+params.UpperMzBand);
      }
      System.out.println("startFloat: "+startFloat);
      System.out.println("stopFloat: "+stopFloat);
      viewer.init();
      viewer.removeSaveLipidomicsSettings();
      
      if (this.displayTolerancePanel_.getShow2D().isSelected()){
        java.awt.Panel l2DPainterPanel = new java.awt.Panel();
        Lipidomics2DPainter l2DPainter = null;
        if (params!=null){
          String[] rawLines2D = rawLines;
          float startFloat2D = startFloat;
          float stopFloat2D = stopFloat;
          if (this.displayTolerancePanel_.getLockMzRange().isSelected() && ((currentIsotopicMass-2*params.LowerMzBand)<startFloat || stopFloat<(currentIsotopicMass+2*params.UpperMzBand))){
            System.out.println("I change the raw lines 2D: "+currentIsotopicMass);
            startFloat2D = currentIsotopicMass-2*params.LowerMzBand;
            stopFloat2D = currentIsotopicMass+2*params.UpperMzBand;
            rawLines2D = reader_.getRawLines(startFloat2D, stopFloat2D, result_.getMsLevels().get(currentSelectedSheet_));
          }
          
          l2DPainter = new Lipidomics2DPainter(analyzer_,rawLines2D, rtTimes, reader_.getRetentionTimes(),
            startFloat2D,stopFloat2D,reader_.getMultiplicationFactorForInt_()/reader_.getLowestResolution_(),params_.LowerMzBand*2,this,
            MSMapViewer.DISPLAY_TIME_MINUTES,currentIsotopicMass-params.LowerMzBand, currentIsotopicMass+params.UpperMzBand, false,
            storedProbes,selectedProbes,Integer.parseInt((String)this.isotope_.getSelectedItem()),charge,result_.getMsLevels().get(currentSelectedSheet_),
            shotgunIsDisplayed_);
          l2DPainter.preChromatogramExtraxtion(currentIsotopicMass-params.LowerMzBand, currentIsotopicMass+params.UpperMzBand);
          l2DPainterPanel = l2DPainter;
        }
        this.makeDisplayRemoveOperations();
        this.viewer_ = viewer;

        majorSplitPane_ = new JSplitPane(JSplitPane.VERTICAL_SPLIT);
        majorSplitPane_.setDividerSize(1);
        majorSplitPane_.setTopComponent(viewer_);
        l2dPanel_.add(l2DPainterPanel,BorderLayout.CENTER);
        majorSplitPane_.setBottomComponent(this.l2dPanel_);
        l2DPainterPanel.setBackground(Color.WHITE);
      
        this.displayPanel_.add(majorSplitPane_,BorderLayout.CENTER);
        this.displayPanel_.validate();
        majorSplitPane_.setDividerLocation(0.75);
      
        if (l2DPainter!=null) {
          l2DPainter_ = l2DPainter;
          l2DPainter_.draw2DDiagram(currentIsotopicMass-params.LowerMzBand, currentIsotopicMass+params.UpperMzBand, m_chkRaw_.isSelected());
        }
        majorSplitPane_.repaint();
        l2DPainterPanel.repaint();
      }else{
        this.makeDisplayRemoveOperations();
        this.viewer_ = viewer;
        this.displayPanel_.add(viewer_,BorderLayout.CENTER);
        this.displayPanel_.validate();
      }
    }
    catch (CgException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    this.showInputElements();
  }
  
  private void makeDisplayRemoveOperations(){
    makeDisplayRemoveOperations(true);
  }
  
  private void makeDisplayRemoveOperations(boolean remove2DPainter){
    this.hideInputElements();
    if (viewer_!=null){
      try{
        viewer_.destroyViewer();
      }catch(Exception ex){
        ex.printStackTrace();
      }
      if (majorSplitPane_!=null)
        majorSplitPane_.remove(viewer_);
      this.displayPanel_.remove(viewer_);
      this.viewer_ = null;
    }
    if (remove2DPainter)
      remove2DPainter();
    if (spectrumPainter_!=null){
      spectrumPainter_.getGraphics().dispose();
      if (finalButtonSection_!=null) spectrumPanel_.remove(finalButtonSection_);
      spectrumPanel_.remove(spectrumPainter_);
      spectrumPainter_ = null;
    }
    if (topSplitPane_!=null){
      topSplitPane_.removeAll();
      topSplitPane_.getGraphics().dispose();
      this.majorSplitPane_.remove(topSplitPane_);
      topSplitPane_ = null;
    }
    if (majorSplitPane_!=null){
      majorSplitPane_.removeAll();
      majorSplitPane_.getGraphics().dispose();
      this.displayPanel_.remove(majorSplitPane_);
      majorSplitPane_= null;
    }
    if(msnUserInterfaceObject_ != null)
      msnUserInterfaceObject_.deleteDetailsBoxes();
    System.gc();
  }


  /**
   * removes graphics methods for displaying the 2D chromatogram/spectrum viewer
   */
  private void remove2DPainter() {
    if (l2DPainter_!=null){
      if (l2DPainter_.getGraphics()!=null)
        l2DPainter_.getGraphics().dispose();
      l2dPanel_.remove(l2DPainter_);
      l2DPainter_ = null;
    }
  }
  
  
  public static void main(String[] args)
  {
    System.out.println("LipidDataAnalyzer "+Settings.VERSION);
    String lookAndFeel = "system";
    try{
      File file = new File(Settings.SETTINGS_FILE);
      FileInputStream inNew = new FileInputStream(file);
      Properties properties = new Properties();
      properties.load(inNew);
      inNew.close();
      lookAndFeel = properties.getProperty("LookAndFeel", "default");
      Settings.isWindows();
      if (lookAndFeel.equalsIgnoreCase("default")){
        if (Settings.isWindows()){
          lookAndFeel = "system";
        }else{
          lookAndFeel = "java";
        }  
      }
    }catch(Exception e){
      e.printStackTrace();
    }
    
    //ToolTipManager.sharedInstance().setLightWeightPopupEnabled(false);
    
    if (lookAndFeel.equalsIgnoreCase("system")){
      try {
        UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
      } catch (Exception unused) {
        ; // Ignore exception because we can't do anything.  Will use default.
      }
    }else if (lookAndFeel.equalsIgnoreCase("java")){
      try {
        UIManager.setLookAndFeel(UIManager.getCrossPlatformLookAndFeelClassName());
      } catch (Exception unused) {
        ; // Ignore exception because we can't do anything.  Will use default.
      }      
    }else if (lookAndFeel.equalsIgnoreCase("windows")){
      try {
        UIManager.setLookAndFeel("com.sun.java.swing.plaf.windows.WindowsLookAndFeel");
      } catch (Exception unused) {
        ; // Ignore exception because we can't do anything.  Will use default.
      }      
    }else if (lookAndFeel.equalsIgnoreCase("motif")){
      try {
        UIManager.setLookAndFeel("com.sun.java.swing.plaf.motif.MotifLookAndFeel");
      } catch (Exception unused) {
        ; // Ignore exception because we can't do anything.  Will use default.
      }      
    }else if (lookAndFeel.equalsIgnoreCase("mac")){
      try {
        UIManager.setLookAndFeel("javax.swing.plaf.mac.MacLookAndFeel");
      } catch (Exception unused) {
        ; // Ignore exception because we can't do anything.  Will use default.
      }      
    }else if (lookAndFeel.equalsIgnoreCase("nimbus")){
      try {
        UIManager.setLookAndFeel("javax.swing.plaf.nimbus.NimbusLookAndFeel");
      } catch (Exception unused) {
        ; // Ignore exception because we can't do anything.  Will use default.
      }      
    }
    
    frame_ = new MainFrame(new LipidDataAnalyzer(), FRAME_WIDTH, FRAME_HEIGHT);
    frame_.setTitle(getFrameTitleString());
  }
  
  private static String getFrameTitleString(){
    String titleString = "Lipid Data Analyzer "+Settings.VERSION+"   "+LipidomicsConstants.getCurrentMSMachine()+" settings ";
    String fragSelected = Settings.getFragmentSettingsString();
    if (fragSelected!=null) titleString += " "+fragSelected;    
    titleString += "         \u00A9 2024 - J\u00fcrgen Hartler, Andreas Ziegl, Gerhard G Thallinger, Leonida M Lamp - GNU GPL v3 license";
    return titleString;
  }

  public void listSelectionChanged(int leadIndex){
    if (currentSelected_!=leadIndex || !currentSelectedSheet_.equalsIgnoreCase((String)selectedSheet_.getSelectedItem())){
      currentSelected_=leadIndex;
      currentSelectedSheet_ = (String)selectedSheet_.getSelectedItem();
      params_ = result_.getIdentifications().get(selectedSheet_.getSelectedItem()).get(resultPositionToOriginalLoopkup_.get(leadIndex));
      if (changeIsotopeListener_!=null) isotope_.removeItemListener(changeIsotopeListener_);
      isotope_.removeAllItems();
      if (params_.getMinIsotope()<0){
        for (int j=0; j!=params_.getMinIsotope()-2; j--){
          isotope_.addItem(String.valueOf(j));
        }
      } else {
        for (int j=0; j!=params_.getMaxIsotope()+2; j++){
          isotope_.addItem(String.valueOf(j));
        }
      }
      isotope_.setSelectedIndex(0);
      isotope_.addItemListener(changeIsotopeListener_);
      this.initMS1OrMS2View(params_);
//      initANewViewer(params_);
    }
  }
    
  private void handleTimerEvent(){
    if (this.rawmzThread_!=null && this.rawmzThread_.finished()){
      if (this.rawmzThread_.getErrorString()!=null&&this.rawmzThread_.getErrorString().length()>0){
        @SuppressWarnings("unused")
        WarningMessage dlg = new WarningMessage(new JFrame(), "Error", rawmzThread_.getErrorString());
        this.singleQuantMenu_.getStartQuantification().setEnabled(true);
        singleQuantMenu_.getSpinnerLabel().setVisible(false);
      }
      this.rawmzThread_ = null;
      this.singleQuantMenu_.getProgressBar().setValue(30);
      this.singleQuantMenu_.getQuantifyingLabel().setText("Translating to chrom");
      System.out.println("Translating to chrom from timer");
      this.readFromRaw_ = true;
      String filePath = singleQuantMenu_.getSelectedMzxmlDirectory().getText();
      String suffix = filePath.substring(filePath.lastIndexOf("."));
      if (suffix.equalsIgnoreCase(".wiff")){
        Vector<File> filesToTranslate = BatchQuantThread.getMzXMLFilesOfWiffConversion(filePath);
        if (filesToTranslate.size()==1){
          singleQuantMenu_.getSelectedMzxmlDirectory().setText(filesToTranslate.get(0).getAbsolutePath());
          mzToChromThread_ = new MzxmlToChromThread(filesToTranslate.get(0).getAbsolutePath(),batchQuantMenu_.getNrProcessorsChrom());
          mzToChromThread_.start();          
        } else {
          Vector<RawQuantificationPairVO> pairs = new Vector<RawQuantificationPairVO>();
          File quantFile = new File(singleQuantMenu_.getSelectedQuantDir().getText());
          for (File rawFile: filesToTranslate){
            pairs.add(new RawQuantificationPairVO(rawFile,quantFile,true));
          }
          int amountOfIsotopes = 0;
          int isotopesMustMatch = 0;
          if (this.singleQuantMenu_.getIsoValidation().isSelected()){
            amountOfIsotopes = Integer.parseInt(this.singleQuantMenu_.getAmountOfIsotopes().getText());
            isotopesMustMatch = Integer.parseInt(this.singleQuantMenu_.getAmountOfMatchingSearchIsotopes().getText());
          }
          boolean ok = true;
          float cutoff = 0f;
          float rtShift = 0f;
          try{
            if (singleQuantMenu_.getCutoff().getText()!=null && singleQuantMenu_.getCutoff().getText().length()>0)
              cutoff = Float.parseFloat(singleQuantMenu_.getCutoff().getText().replaceAll(",", "."));
          }catch(NumberFormatException ex){new WarningMessage(new JFrame(), "Error", "The cutoff value must be float format!"); ok=false;}
          try{
            if (singleQuantMenu_.getRtShift().getText()
            		!=null && singleQuantMenu_.getRtShift().getText().length()>0)
              rtShift = Float.parseFloat(singleQuantMenu_.getRtShift().getText().replaceAll(",", "."));
          }catch(NumberFormatException ex){new WarningMessage(new JFrame(), "Error", "The RT-shift value must be float format!"); ok=false;}
          if (ok){
            boolean ionMode=false;
            if (this.singleQuantMenu_.getIonMode()!=null && ((String)singleQuantMenu_.getIonMode().getSelectedItem()).equalsIgnoreCase("+"))
              ionMode = true;
            batchQuantMenu_.getBatchQuantTableModel().clearFiles();
            batchQuantMenu_.getBatchQuantTableModel().addFiles(pairs);
            batchQuantThread_ = new BatchQuantThread(singleQuantMenu_,
                amountOfIsotopes,isotopesMustMatch, cutoff, 
                rtShift,ionMode,false);
            this.batchQuantMenu_.getQuantifyingLabel().setText("Quantifying");
            this.batchQuantMenu_.getProgressBar().setValue(0);
            this.batchQuantMenu_.getQuantifyingPanel().setVisible(true);
            this.batchQuantMenu_.getStartQuantification().setEnabled(false);
            this.batchQuantMenu_.getSpinnerLabel().setVisible(true);
            
            quantitationPane_.setSelectedIndex(0);
            batchQuantThread_.start();
          }
        }
      }else{
        String mzXMLFilePath = singleQuantMenu_.getSelectedMzxmlDirectory().getText().substring(0,singleQuantMenu_.getSelectedMzxmlDirectory().getText().lastIndexOf("."))+"."+LipidomicsConstants.getIntermediateFileFormat();
        mzToChromThread_ = new MzxmlToChromThread(mzXMLFilePath,batchQuantMenu_.getNrProcessorsChrom());
        mzToChromThread_.start();
      }
    }
    if (this.mzToChromThread_!=null && this.mzToChromThread_.finished()){
      if (this.mzToChromThread_.getErrorString()!=null&&this.mzToChromThread_.getErrorString().length()>0){
        @SuppressWarnings("unused")
        WarningMessage dlg = new WarningMessage(new JFrame(), "Error", mzToChromThread_.getErrorString());
        this.singleQuantMenu_.getStartQuantification().setEnabled(true);
        singleQuantMenu_.getSpinnerLabel().setVisible(false);
      }
      if (readFromRaw_){
        this.readFromRaw_ = false;
        RawToMzxmlThread.deleteMzXMLFiles(singleQuantMenu_.getSelectedMzxmlDirectory().getText().substring(0,singleQuantMenu_.getSelectedMzxmlDirectory().getText().lastIndexOf("."))+"."+LipidomicsConstants.getIntermediateFileFormat());
      }

      boolean ok = true;
      if (mzToChromThread_.isPolaritySwitched()){
        if (singleQuantMenu_.getSelectedQuantDir().getText().indexOf(GlobalConstants.CHROMATOGRAM_HEADER_FILE_POLARITY_POSITIVE)!=-1){
          String newRawName = singleQuantMenu_.getSelectedMzxmlDirectory().getText();
          newRawName = newRawName.substring(0,newRawName.lastIndexOf("."))+RawToChromThread.FILE_SUFFIX_POLARITY_POSITIVE+".chrom";
          singleQuantMenu_.getSelectedMzxmlDirectory().setText(newRawName);
        } else if (singleQuantMenu_.getSelectedQuantDir().getText().indexOf(GlobalConstants.CHROMATOGRAM_HEADER_FILE_POLARITY_NEGATIVE)!=-1){
          String newRawName = singleQuantMenu_.getSelectedMzxmlDirectory().getText();
          newRawName = newRawName.substring(0,newRawName.lastIndexOf("."))+RawToChromThread.FILE_SUFFIX_POLARITY_NEGATIVE+".chrom";
          singleQuantMenu_.getSelectedMzxmlDirectory().setText(newRawName);
        } else {
          ok = false;
          new WarningMessage(new JFrame(), "Error", "This is polarity switched data! It is not clear which of the two generated chrom files shall be quantified! Please use "+GlobalConstants.CHROMATOGRAM_HEADER_FILE_POLARITY_POSITIVE+" or "+GlobalConstants.CHROMATOGRAM_HEADER_FILE_POLARITY_NEGATIVE+" in the file name!");
          this.singleQuantMenu_.getQuantifyingPanel().setVisible(false);
        }
      }
      this.mzToChromThread_ = null;
      this.singleQuantMenu_.getQuantifyingLabel().setText("Quantifying");
      System.out.println("Quantifying from thread");
      this.singleQuantMenu_.getProgressBar().setValue(70);
      int amountOfIsotopes = 0;
      int isotopesMustMatch = 0;
      if (this.singleQuantMenu_.getIsoValidation().isSelected()){
        amountOfIsotopes = Integer.parseInt(this.singleQuantMenu_.getAmountOfIsotopes().getText());
        isotopesMustMatch = Integer.parseInt(this.singleQuantMenu_.getAmountOfMatchingSearchIsotopes().getText());
      }      
      float cutoff = 0f;
      float rtShift = 0f;
      try{
        if (singleQuantMenu_.getCutoff().getText()!=null && singleQuantMenu_.getCutoff().getText().length()>0)
          cutoff = Float.parseFloat(singleQuantMenu_.getCutoff().getText().replaceAll(",", "."));
      }catch(NumberFormatException ex){new WarningMessage(new JFrame(), "Error", "The cutoff value must be float format!"); ok=false;}
      try{
        if (singleQuantMenu_.getRtShift().getText()!=null && singleQuantMenu_.getRtShift().getText().length()>0)
          rtShift = Float.parseFloat(singleQuantMenu_.getRtShift().getText().replaceAll(",", "."));
      }catch(NumberFormatException ex){new WarningMessage(new JFrame(), "Error", "The RT-shift value must be float format!"); ok=false;}
      if (ok){
        boolean ionMode = false;
        if (this.singleQuantMenu_.getIonMode()!=null && ((String)singleQuantMenu_.getIonMode().getSelectedItem()).equalsIgnoreCase("+"))
          ionMode = true;
        quantThread_ = new QuantificationThread(singleQuantMenu_.getSelectedMzxmlDirectory().getText(), singleQuantMenu_.getSelectedQuantDir().getText(),LipidDataAnalyzer.getResultFilePath(singleQuantMenu_.getSelectedMzxmlDirectory().getText(), singleQuantMenu_.getSelectedQuantDir().getText()),// Float.parseFloat(this.singleMzTol_.getText()),
            Float.parseFloat(this.singleQuantMenu_.getTimeMinusTol().getText()),Float.parseFloat(this.singleQuantMenu_.getTimePlusTol().getText()),amountOfIsotopes,isotopesMustMatch,this.singleQuantMenu_.getSearchUnknownTime().isSelected(),
            cutoff,rtShift,this.singleQuantMenu_.getNrProcessors(),ionMode,false);
        quantThread_.start();
      }else{
        this.singleQuantMenu_.getStartQuantification().setEnabled(true);
        singleQuantMenu_.getSpinnerLabel().setVisible(false);
      }
    }
    if (this.quantThread_!=null && !this.quantThread_.finished() && 
        (this.quantThread_.getErrorString()==null||this.quantThread_.getErrorString().length()==0)){
      if (quantThread_.getTotalAmountOfLipids()>0&&quantThread_.getCurrentLipidCount()>0){
        this.singleQuantMenu_.getProgressBar().setValue(70+((30*(quantThread_.getCurrentLipidCount()-1))/quantThread_.getTotalAmountOfLipids()));
        singleQuantMenu_.getQuantifyingLabel().setText("Quantifying "+quantThread_.getCurrentLipid()+" ("+quantThread_.getCurrentLipidCount()+"/"+quantThread_.getTotalAmountOfLipids()+")");
      }  
    }
    if (this.quantThread_!=null && this.quantThread_.finished()){
      if (this.quantThread_.getErrorString()!=null&&this.quantThread_.getErrorString().length()>0){
        @SuppressWarnings("unused")
        WarningMessage dlg = new WarningMessage(new JFrame(), "Error", quantThread_.getErrorString());
        this.singleQuantMenu_.getStartQuantification().setEnabled(true);
        singleQuantMenu_.getSpinnerLabel().setVisible(false);
      }
      this.quantThread_ = null;
      this.singleQuantMenu_.getProgressBar().setValue(100);
      singleQuantMenu_.getQuantifyingLabel().setText("Finished");
      this.singleQuantMenu_.getStartQuantification().setEnabled(true);
      singleQuantMenu_.getSpinnerLabel().setVisible(false);
    }
    if (this.batchQuantThread_!=null && this.batchQuantThread_.finished()){
      this.batchQuantThread_ = null;
      this.batchQuantMenu_.getProgressBar().setValue(100);
      batchQuantMenu_.getQuantifyingLabel().setText("Finished");
      this.batchQuantMenu_.getStartQuantification().setEnabled(true);
      batchQuantMenu_.getSpinnerLabel().setVisible(false);
      
      if (this.singleQuantMenu_.getProgressBar().getValue() > 0) //possible with .wiff files
      {
      	this.singleQuantMenu_.getProgressBar().setValue(100);
        singleQuantMenu_.getQuantifyingLabel().setText("Finished");
        this.singleQuantMenu_.getStartQuantification().setEnabled(true);
        singleQuantMenu_.getSpinnerLabel().setVisible(false);
      }
    }

  }
  
  public static String getResultFilePath(String rawFilePath, String quantFilePath){
    String resultFilePath = rawFilePath.substring(0,rawFilePath.lastIndexOf("."));
    String quantFileFrag = new String(quantFilePath);
    int index = -1;
    int slashIndex = quantFileFrag.lastIndexOf("/");
    int backSlashIndex = quantFileFrag.lastIndexOf("\\");
    if (slashIndex>backSlashIndex)
      index = slashIndex;
    else
      index = backSlashIndex;
    resultFilePath+="_"+quantFilePath.substring(index+1);
    if (resultFilePath.indexOf(".")!=-1) resultFilePath=resultFilePath.substring(0,resultFilePath.lastIndexOf("."));
    resultFilePath += ".xlsx";
    return resultFilePath;
  }
  
  private class ThreadSupervisor extends TimerTask{

    public void run()
    {
      handleTimerEvent();
    }
    
  }
  
  private void changeSeletectSheet(String command){
    this.updateResultListSelectionTable();
  }
  
  private void updateResultListSelectionTable(){
    tablePanel_.remove(tablePane);
    
    int orderType = LipidomicsJTable.ORDER_TYPE_AS_IS;
    if (orderResultsType_.containsKey(selectedSheet_.getSelectedItem()))
      orderType = orderResultsType_.get(selectedSheet_.getSelectedItem());
    resultPositionToOriginalLoopkup_ = new Hashtable<Integer,Integer>();
    resultPositionToMolecularSpeciesLookup_ = new Hashtable<Integer,String>();
    Vector<LipidParameterSet> lipids = result_.getIdentifications().get(selectedSheet_.getSelectedItem());
    Vector<LipidParameterSet> lipidsOrdered = new Vector<LipidParameterSet>();
    if (orderType==LipidomicsJTable.ORDER_TYPE_MZ || orderType==LipidomicsJTable.ORDER_TYPE_INTENSITY){
      for (LipidParameterSet analyte : lipids){
        int count = 0;
        for (int i=0;i!=lipidsOrdered.size();i++){
          if (orderType==LipidomicsJTable.ORDER_TYPE_MZ){
            if (lipidsOrdered.get(i).Mz[0]>analyte.Mz[0]) break;
          }else if (orderType==LipidomicsJTable.ORDER_TYPE_INTENSITY){
            if (analyte.getArea()>lipidsOrdered.get(i).getArea()) break;
          }
          count++;
        }
        lipidsOrdered.add(count,analyte);
      }
    } else lipidsOrdered = new Vector<LipidParameterSet>(lipids);
    LipidomicsTableModel model = new LipidomicsTableModel(lipidsOrdered,lipids,displayTolerancePanel_.getShowMSnNames().isSelected(),displayTolerancePanel_.getShowOmegaNames().isSelected(),resultsShowModification_.get(selectedSheet_.getSelectedItem()));
    resultPositionToOriginalLoopkup_ = model.getPositionToOriginal(); 
    resultPositionToMolecularSpeciesLookup_ = model.getRowToMolecularSpeciesHumanReadable();
    displayTable_ = new LipidomicsJTable(model, new LipidomicsTableCellRenderer(),
        reader_.getHighestMsLevel()>1&&reader_.getMsmsType().equalsIgnoreCase(ChromatogramReader.CHROMATOGRAM_HEADER_FILE_MSMS_TYPE_PRECURSOR), orderType,
        QuantificationThread.hasRtInfo(result_.getIdentifications()),this);
    listSelectionModel = displayTable_.getSelectionModel();
    displayTable_.setSelectionModel(listSelectionModel);
    listSelectionModel.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
    tablePane = new JScrollPane(displayTable_);
    tablePane.setPreferredSize(new Dimension(420, 130));
    tablePanel_.add(tablePane,BorderLayout.CENTER);
    tablePanel_.invalidate();
    tablePanel_.updateUI();
    selectionPane.setVisible(true);    
  }
  
  private void updateAnalysisSelectionTable(){
    analysisSelectionTablePanel_.remove(analysisTablePane_);
    String[][] tableData = new String[resultFiles_.size()][2];
    int count=0;
    for (File file : resultFiles_){
      String fileName = StaticUtils.extractFileName(file.getAbsolutePath()); 
      String dir = StaticUtils.extractDirName(file.getAbsolutePath()); 
      tableData[count][0] = fileName;
      tableData[count][1] = dir;
      count++;
    }
    this.generateResultsAnalysisTablePane(tableData);
    analysisSelectionTablePanel_.invalidate();
    analysisSelectionTablePanel_.updateUI();
  }
  
  private void generateResultsAnalysisTablePane(String[][] tableData){
    int tableWidth = RESULTS_TABLE_WIDTH;
    //18 is the width of the scrollbar and needs to be subtracted from tableWidth to get the space available for the columns
    int columnWidth = tableWidth-18;
    String[] columnNames = { "File Name", "Directory" };
    resultFilesDisplayTable_ = new JTable(tableData, columnNames);
    resultListSelectionModel_ = resultFilesDisplayTable_.getSelectionModel();
    resultFilesDisplayTable_.setSelectionModel(resultListSelectionModel_);
    resultListSelectionModel_.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);  
    analysisTablePane_ = new JScrollPane(resultFilesDisplayTable_);
    analysisTablePane_.setPreferredSize(new Dimension(tableWidth, 300));
    if (tableData!=null && tableData.length>0){
      int longestFirstColumnElement = 0;
      int longestSecondColumnElement = 0;
      for (int i=0;i!=tableData.length;i++){
        int lengthOne = tableData[i][0].length();
        int lengthTwo = tableData[i][1].length();
        if (lengthOne>longestFirstColumnElement)longestFirstColumnElement = lengthOne;
        if (lengthTwo>longestSecondColumnElement)longestSecondColumnElement = lengthTwo;
      }
      double percentOne = (double)longestFirstColumnElement/(double)(longestFirstColumnElement+longestSecondColumnElement);
      if (percentOne<0.25)
        percentOne = 0.25;
      if (percentOne>0.75)
        percentOne = 0.75;
      int columnWidthOne = (int)((double)columnWidth*percentOne);
      int columnWidthTwo = columnWidth-columnWidthOne;
      resultFilesDisplayTable_.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
      resultFilesDisplayTable_.getColumnModel().getColumn(0).setPreferredWidth(columnWidthOne);
      resultFilesDisplayTable_.getColumnModel().getColumn(1).setPreferredWidth(columnWidthTwo);
    }
    analysisSelectionTablePanel_.add(analysisTablePane_);
    
  }
   
  
  protected void change2DTypeState(ItemEvent e, String command){
    if (command.equalsIgnoreCase("DisplayModeRaw") && this.l2DPainter_!=null){
      this.l2DPainter_.setRaw(true);
      this.l2DPainter_.repaint();
    }
    if (command.equalsIgnoreCase("DisplayModeSmooth") && this.l2DPainter_!=null){
      this.l2DPainter_.setRaw(false);
      this.l2DPainter_.repaint();
    }
    if (command.equalsIgnoreCase(DisplayTolerancePanel.SHOW_2D_CHANGED)){
      if (this.params_!=null){
        if (this.displayTolerancePanel_.getShow2D().isSelected())
          this.initMS1OrMS2View(params_);
        else
          this.initANewViewer(params_);
      }
///        this.initANewViewer(params_);
    }
    if (command.equalsIgnoreCase(DisplayTolerancePanel.SHOW_MSN_ONLY)){
    	if(this.displayTolerancePanel_.getShowMSnEvidence().isSelected())
    	{
    		result_ = mSnResult_;
    		displayTolerancePanel_.getShowChainEvidence().setEnabled(true);
    	}
    	else {
    		displayTolerancePanel_.getShowChainEvidence().setSelected(false);
    		displayTolerancePanel_.getShowChainEvidence().setEnabled(false);
    		result_ = originalResult_;
    	}
      this.updateResultListSelectionTable();
    }
    if (command.equalsIgnoreCase(DisplayTolerancePanel.SHOW_CHAIN_ONLY)){
    	if(this.displayTolerancePanel_.getShowChainEvidence().isSelected())
    	{
    		result_ = chainResult_;
    	}
    	else {
    		result_ = mSnResult_;
    	}
      this.updateResultListSelectionTable();
    }
    if (command.equalsIgnoreCase("showMSnOnlyStat")){
    	if(this.showMSnEvidenceStat_.isSelected())
    	{
    		statisticsViewMode_ = 1;
    		showChainEvidenceStat_.setEnabled(true);
    	}
    	else {
    		showChainEvidenceStat_.setSelected(false);
    		statisticsViewMode_ = 0; //this must be after deselecting showing chain evidence, as this has a side effect on the statistics view mode.
    		showChainEvidenceStat_.setEnabled(false);
    	}
    }
    if (command.equalsIgnoreCase("showChainOnlyStat")){
    	if(this.showChainEvidenceStat_.isSelected())
    	{
    		statisticsViewMode_ = 2;
    	}
    	else {
    		statisticsViewMode_ = 1; //if this option is not selected, the default is that MSn only is selected.
    	}
    }
    if (command.equalsIgnoreCase("combineWithOx")){
    	if(this.combineWithOx_.isSelected())
    	{
    		combineOxWithNonOx_ = true;
    	}
    	else {
    		combineOxWithNonOx_ = false;
    	}
    }
    if (command.equalsIgnoreCase("DisplayModeAbundance")){
      this.spectrumPainter_.setRelativeValues(relAbund_.isSelected());
      this.spectrumPainter_.repaint();
    }
    if (command.equalsIgnoreCase("endisableAmountOfIsotopes")){
      if (this.singleQuantMenu_.getIsoValidation().isSelected()){
        this.singleQuantMenu_.getAmountOfIsotopes().setEnabled(true);
        this.singleQuantMenu_.getAmountOfMatchingSearchIsotopes().setEnabled(true);
      }else{
        this.singleQuantMenu_.getAmountOfIsotopes().setEnabled(false);
        this.singleQuantMenu_.getAmountOfMatchingSearchIsotopes().setEnabled(false);
      }  
    }
    if (command.equalsIgnoreCase("endisableAmountOfBatchIsotopes")){
      if (this.batchQuantMenu_.getIsoValidation().isSelected()){
        this.batchQuantMenu_.getAmountOfIsotopes().setEnabled(true);
        this.batchQuantMenu_.getAmountOfMatchingSearchIsotopes().setEnabled(true);
      }else{
        this.batchQuantMenu_.getAmountOfIsotopes().setEnabled(false);
        this.batchQuantMenu_.getAmountOfMatchingSearchIsotopes().setEnabled(false);
      }  
    }
//    if (command.equalsIgnoreCase("searchForUnknownRetentionTime")){
//      if (this.singleQuantMenu_.getSearchUnknownTime().isSelected())
//        this.singleQuantMenu_.getAmountOfMatchingSearchIsotopes().setEnabled(true);
//      else
//        this.singleQuantMenu_.getAmountOfMatchingSearchIsotopes().setEnabled(false);
//    }
//    if (command.equalsIgnoreCase("searchForUnknownBatchRetentionTime")){
//      if (this.batchQuantMenu_.getSearchUnknownTime().isSelected())
//        this.batchQuantMenu_.getAmountOfMatchingSearchIsotopes().setEnabled(true);
//      else
//        this.batchQuantMenu_.getAmountOfMatchingSearchIsotopes().setEnabled(false);
//    }   
    if (command.equalsIgnoreCase("ChangeIsotope") && !this.displaysMs2_){
      if (e.getStateChange()==ItemEvent.SELECTED){
        if (this.params_!=null && l2DPainter_!=null){
          if (!this.displayTolerancePanel_.getLockMzRange().isSelected() || this.viewer_==null || lockRangeUpdateRequired_){
            this.initANewViewer(params_,this.l2DPainter_.getAllSelectedProbes());
            lockRangeUpdateRequired_ = false;
          }else {
            updateViewForLockMz(params_);
          }        
        }
      }
    } else if (command.equalsIgnoreCase(DisplayTolerancePanel.SHOW_MSN_NAMES)){
      updateResultListSelectionTable();
    } else if (command.equalsIgnoreCase(DisplayTolerancePanel.SHOW_OMEGA_NAMES)) {
      updateResultListSelectionTable();
    } else if (command.equalsIgnoreCase(DisplayTolerancePanel.CHANGE_LOCK_MZ_RANGE)) {
      changeLockMzRange();
    }
    
  }
  
  protected class LipidomicsItemListener implements java.awt.event.ItemListener
  {
    private String m_ctrl;
    
    public LipidomicsItemListener(String ctrl){
      m_ctrl = ctrl;
    }

    public void itemStateChanged(ItemEvent e)  {
      change2DTypeState(e,m_ctrl);
    }
  }
  
  private class SheetSelectionItemListener implements java.awt.event.ItemListener
  {
    private String m_ctrl;
    
    public SheetSelectionItemListener(String ctrl){
      m_ctrl = ctrl;
    }

    public void itemStateChanged(ItemEvent e)  {
      if (m_ctrl.equalsIgnoreCase("ChangeSheet")){
        if (e.getStateChange()==ItemEvent.SELECTED){
          changeSeletectSheet(m_ctrl);
        }
      }
    }
  }
  
  private class FragSettingsChangeListener implements java.awt.event.ItemListener
  {
    private String m_ctrl;
    
    public FragSettingsChangeListener(String ctrl){
      m_ctrl = ctrl;
    }

    public void itemStateChanged(ItemEvent e)  {
      if (m_ctrl.equalsIgnoreCase("ChangeMachine")){
        if (e.getStateChange()==ItemEvent.SELECTED){
          refreshFragSettingsSelection();
        }
      }
      if (m_ctrl.equalsIgnoreCase("ChangeFragSelection")){
        if (e.getStateChange()==ItemEvent.SELECTED){
          refreshFragSettingsVisibility();
        }
      }
    }
  }

  private void refreshFragSettingsSelection(){
    String currentMachine = (String)msMachineTypes_.getSelectedItem();
    fragmentationSettings1_.removeAllItems();
    fragmentationSettings2_.removeAllItems();
    for (String fragSetting : Settings.getFragmentationSettings(currentMachine)) fragmentationSettings1_.addItem(fragSetting);
    for (String fragSetting : Settings.getFragmentationSettings(currentMachine)) fragmentationSettings2_.addItem(fragSetting);
    if (currentMachine.equalsIgnoreCase(LipidomicsConstants.getCurrentMSMachine())){
      fragmentationSettings1_.setSelectedItem(Settings.getFragmentationSelection1());
      fragmentationSettings2_.setSelectedItem(Settings.getFragmentationSelection2());
    } else {
      fragmentationSettings1_.setSelectedItem(Settings.FRAG_SELECTION_NONE);
      fragmentationSettings2_.setSelectedItem(Settings.FRAG_SELECTION_NONE);
    }
    refreshFragSettingsVisibility();
  }
  
  private void refreshFragSettingsVisibility(){
    String selected1 = (String)fragmentationSettings1_.getSelectedItem();
    if (selected1.equalsIgnoreCase(Settings.FRAG_SELECTION_NONE)||selected1.equalsIgnoreCase(Settings.FRAG_SELECTION_NO_INTENSITY))
      inputFragSettings2_.setVisible(false);
    else
      inputFragSettings2_.setVisible(true);
    settingsPanel_.invalidate();
  }
  
  private void hideInputElements(){
    isotopeLabel_.setVisible(false);
    m_chkRaw_.setVisible(false);
    m_chkSmooth_.setVisible(false);
    m_upButton_.setVisible(false);
    m_dnButton_.setVisible(false);
    lx_min_.setVisible(false);
    lx_max_.setVisible(false);
    m_minTimeText_.setVisible(false);
    m_maxTimeText_.setVisible(false);
    spectrumSelectedLabel_.setVisible(false);
    spectrumSelected_.setVisible(false);
    rtSelectedLabel_.setVisible(false);
    rtSelected_.setVisible(false);
    precursorSelectedLabel_.setVisible(false);
    for (JLabel label : precursorSelected_) label.setVisible(false);
    msLevelSelectedLabel_.setVisible(false);
    msLevelSelected_.setVisible(false);
    spectrumEarlier_.setVisible(false);
    spectrumLater_.setVisible(false);
    lz_min_.setVisible(false);
    lz_max_.setVisible(false);
    mz_minTimeText_.setVisible(false);
    mz_maxTimeText_.setVisible(false);
    m_zoomIn_.setVisible(false);
    m_zoomAll_.setVisible(false);
    mz_zoomIn_.setVisible(false);
    mz_zoomAll_.setVisible(false);
    storeSelectedAreas_.setVisible(false);
    isotope_.setVisible(false);
    relAbund_.setVisible(false);
    absAbund_.setVisible(false);
    annotationLabel_.setVisible(false);
    annotationThreshold_.setVisible(false);
    annotationUnit_.setVisible(false);
    if (!exportChromatogramsFromDRView_)
    	exportSpectra_.setVisible(false);
  }
  
  private void showInputElements(){
    isotopeLabel_.setVisible(true);
    m_chkRaw_.setVisible(true);
    m_chkSmooth_.setVisible(true);
    m_upButton_.setVisible(true);
    m_dnButton_.setVisible(true);
    lx_min_.setVisible(true);
    lx_max_.setVisible(true);
    m_minTimeText_.setVisible(true);
    m_maxTimeText_.setVisible(true);
    m_zoomIn_.setVisible(true);
    m_zoomAll_.setVisible(true);
    storeSelectedAreas_.setVisible(true);
    isotope_.setVisible(true);
  }
  
  private void showMs2InputElements(){
    spectrumSelectedLabel_.setVisible(true);
    spectrumSelected_.setVisible(true);
    spectrumEarlier_.setVisible(true);
    spectrumLater_.setVisible(true);
    rtSelectedLabel_.setVisible(true);
    rtSelected_.setVisible(true);
    precursorSelectedLabel_.setVisible(true);
    for (JLabel label : precursorSelected_) label.setVisible(true);
    msLevelSelectedLabel_.setVisible(true);
    msLevelSelected_.setVisible(true);
    relAbund_.setVisible(true);
    absAbund_.setVisible(true);
    lz_min_.setVisible(true);
    lz_max_.setVisible(true);
    mz_minTimeText_.setVisible(true);
    mz_maxTimeText_.setVisible(true);
    mz_zoomIn_.setVisible(true);
    mz_zoomAll_.setVisible(true);
    annotationLabel_.setVisible(true);
    annotationThreshold_.setVisible(true);
    annotationUnit_.setVisible(true);
    showChromInputElementsForMs2();
  }
  
  private void showChromInputElementsForMs2(){
    m_chkRaw_.setVisible(true);
    m_chkSmooth_.setVisible(true);
    lx_min_.setVisible(true);
    lx_max_.setVisible(true);
    m_minTimeText_.setVisible(true);
    m_maxTimeText_.setVisible(true);
    m_zoomIn_.setVisible(true);
    m_zoomAll_.setVisible(true);
    exportSpectra_.setVisible(true);
  }


  private void storeSelectedAreasToFile(){
    Vector<CgProbe> probesToStore = new Vector<CgProbe>();
    Vector<Vector<CgProbe>> allProbes = this.l2DPainter_.getAllSelectedProbes();
    probesToStore.addAll(allProbes.get(0));
    probesToStore.addAll(allProbes.get(1));
    float totalArea = 0;
    int highestIsotope = 0;
    int lowestIsotope = 0;
    for (CgProbe aProbe: probesToStore){
      if (aProbe.AreaStatus==CgAreaStatus.OK){
        totalArea+=aProbe.Area;
        if (aProbe.isotopeNumber>highestIsotope)
          highestIsotope = aProbe.isotopeNumber;
        if (aProbe.isotopeNumber<lowestIsotope)
          lowestIsotope = aProbe.isotopeNumber;
          
      }
    }
    Hashtable<Integer,Boolean> foundIsotopes = new Hashtable<Integer,Boolean>();
    for (int i=lowestIsotope;i!=highestIsotope;i++){
      foundIsotopes.put(i, false);
    }
    for (CgProbe aProbe: probesToStore){
      if (aProbe.AreaStatus==CgAreaStatus.OK && foundIsotopes.containsKey(aProbe.isotopeNumber))
        foundIsotopes.put(aProbe.isotopeNumber, true);
    }
    List<Integer> keyList = new ArrayList<Integer>(foundIsotopes.keySet());
    Collections.sort(keyList);
    Vector<Integer> missingIsotopes = new Vector<Integer>();
    for (Integer isoKey : keyList){
      if (!foundIsotopes.get(isoKey)){
        if (lowestIsotope<0)
          missingIsotopes.add(0,isoKey);
        else
          missingIsotopes.add(isoKey);
      }  
    }
    if (missingIsotopes.size()==0){
      this.result_.getIdentifications().get(this.currentSelectedSheet_).get(resultPositionToOriginalLoopkup_.get(this.currentSelected_)).setProbes(probesToStore);
      this.result_.getIdentifications().get(this.currentSelectedSheet_).get(resultPositionToOriginalLoopkup_.get(this.currentSelected_)).Area = totalArea;
      this.params_ = this.result_.getIdentifications().get(this.currentSelectedSheet_).get(resultPositionToOriginalLoopkup_.get(this.currentSelected_));
      
      try{
        QuantificationResultExporter.writeResultsToExcel(selectedResultFile.getText(), result_);
        this.readResultFile(selectedResultFile.getText(),true);
        this.updateResultListSelectionTable();
        this.displayTable_.changeSelection(this.currentSelected_, 1, false, false);
        listSelectionChanged(this.currentSelected_);
        this.l2DPainter_.setStoredProbes(probesToStore);
        allProbes = new Vector<Vector<CgProbe>>();
        allProbes.add(probesToStore);
        allProbes.add(new Vector<CgProbe>());
        this.viewer_.setPaintableProbes(allProbes);
        this.viewer_.setTheShowSelectedWasOn(false);
        this.viewer_.repaintColors();
      }catch (Exception ex){
        ex.printStackTrace();
        @SuppressWarnings("unused")
        WarningMessage dlg = new WarningMessage(new JFrame(), "Error", "The areas could not be stored: "+ex.getMessage());
      }
    }else{
      String missingString = "";
      for (Integer missIso : missingIsotopes){
        missingString+=" ,"+missIso;
      }
      missingString = missingString.substring(2);
      new WarningMessage(new JFrame(), "Error", "You cannot store a +"+highestIsotope+" isotope without having selected peaks at "+missingString+" isotope! So add these isotopes or delete the +"+highestIsotope+" isotope!");
    }
  }
  
  private void updateSheetSelectionList(){
    
    tablePanel_.remove(selectedSheet_);
    selectedSheet_ = new JComboBox<String>();
    selectedSheet_.addItemListener(new SheetSelectionItemListener("ChangeSheet"));
    selectedSheet_.setPreferredSize(new Dimension(30, 20));
    
    //Sort the keys for oxLipids
    List<String> sortedSheets = Collections.list(result_.getIdentifications().keys());
    Collections.sort(sortedSheets);
    
    for (String sheetName : sortedSheets) {
      selectedSheet_.addItem(sheetName);
    }
    
//    for (String sheetName : result_.getIdentifications().keySet()){ 
//      selectedSheet_.addItem(sheetName);
//    }
    selectedSheet_.setToolTipText(TooltipTexts.DISPLAY_SELECT_CLASS);
    tablePanel_.add(selectedSheet_,BorderLayout.NORTH);
    tablePanel_.invalidate();
    tablePanel_.updateUI();
    selectionPane.setVisible(true);
  }

  public boolean heatMapClicked(String experimentName, ResultCompVO compVO, String moleculeName, boolean showMSn)
  {
    File resultsFile = new File(compVO.getAbsoluteFilePath());
    if (resultsFile.exists()&&resultsFile.isFile()){
      selectedResultFile.setText(compVO.getAbsoluteFilePath());
      String chromFileBase = StaticUtils.extractChromBaseName(compVO.getAbsoluteFilePath(),experimentName);
    	try
    	{
    		this.readResultFile(selectedResultFile.getText());
    	}
    	catch (ExcelInputFileException ex){}
    	
      if (chromFileBase == null)
      {
      	String rawName = result_.getConstants().getRawFileName();
      	String filePath = resultsFile.getParentFile().toString();
      	chromFileBase = filePath+"\\"+rawName.substring(0,rawName.indexOf("."));
      }
      boolean chromFileExists = false;
      if (chromFileBase!=null && chromFileBase.length()>0){
        chromFileExists = true;        
      }
      if (chromFileExists){
        if (!StaticUtils.existChromFiles(chromFileBase)){
          return false;
        }
        String chromFilePath = chromFileBase+".chrom";
        System.out.println("chromFilePath: "+chromFilePath);
        this.selectedChromFile.setText(chromFilePath);      
        startDisplay();
        String sheetToSelect = resultTabs_.getTitleAt(resultTabs_.getSelectedIndex());
        selectedSheet_.setSelectedItem(sheetToSelect);
        displayTolerancePanel_.getShowMSnNames().setSelected(showMSn);
        displayTolerancePanel_.getShowOmegaNames().setSelected(displayTolerancePanel_.isShowOmegaEditingTools() && displayTolerancePanel_.getShowMSnNames().isSelected());
        String moleculeInTableName = null;
        int selection = -1;
        for (int i=0;i!=this.displayTable_.getRowCount();i++){
          String moleculeInTable = (String)this.displayTable_.getDisplayedNameAt(i);
          String rtInTableString = null;
          if (moleculeInTable.startsWith(moleculeName)){
            boolean found = false;
            if (moleculeInTable.equals(moleculeName)) {
            	found = true;
            }else {
            	rtInTableString = moleculeInTable.substring(moleculeName.length()+1);
            	if (rtInTableString.indexOf("_")!=-1) rtInTableString = rtInTableString.substring(0,rtInTableString.indexOf("_"));
            	if (compVO.getResultMolecule().belongsRtToThisAreaVO(rtInTableString, null))
            	{
            		found = true;
            	}
            }
            if (found){
              //if show MS2 spectra is selected, try to find an adequate matching hit where MS2 spectra are present
              if (this.displaysMs2_){
                int j=i;
                boolean foundMsn = false;
                while (j<displayTable_.getRowCount() && ((String)this.displayTable_.getDisplayedNameAt(j)).startsWith(moleculeName)){
                  if (((LipidomicsTableModel)displayTable_.getModel()).hasMS2Evidence(j)){
                    if (compVO.getResultMolecule().belongsRtToThisAreaVO(rtInTableString, null))
                    {
                    	foundMsn = true;
                    	break;
                    }
                  }
                  j++;
                }
                if (foundMsn && i!=j){
                  i = j;
                  moleculeInTable = (String)this.displayTable_.getDisplayedNameAt(i);
                }
              }
              moleculeInTableName = moleculeInTable;
              selection = i;
              break;
            }
          }
        }
        if (moleculeInTableName!=null){
        	mainTabs.setSelectedIndex(mainTabs.indexOfComponent(displayPanel_));
//          try {
//            Thread.sleep(100);
//          }
//          catch (InterruptedException e) {
//            e.printStackTrace();
//          }
//          ListSelectionEvent event2 = new ListSelectionEvent(displayTable,selection,selection+1,false);
//          displayTable.addSelectionInterval(selection,selection);
//          System.out.println("Handlers: "+handler.length);
//          .valueChanged(event2);
          this.makeDisplayRemoveOperations();
          this.displayTable_.changeSelection(selection, 1, false, false);
          listSelectionChanged(selection);
        }else
          return false;
      }else{
        new WarningMessage(new JFrame(),"ERROR","The chrom file \""+experimentName+".chrom\" is not there!");
        return false;
      }  
    }else{
      new WarningMessage(new JFrame(),"ERROR","The result file \""+compVO.getAbsoluteFilePath()+"\" does not exist!");
      return false;
    }  
    return true;
  }

  
  public boolean analyteClicked(String moleculeName, String groupName, int maxIsotope, boolean rtGrouped, ResultDisplaySettingsVO settingVO, String prefUnit, String unit){
    molBarCharts_.get(groupName).remove(1);
    Hashtable<String,ResultCompVO> analysisResults = new Hashtable<String,ResultCompVO>(analysisModule_.getResults().get(groupName).get(moleculeName));
//    Vector<String> expNames = new Vector<String>();
//    for (String name : analysisModule_.getExpNamesInSequence()){
//      String displayName  = getDisplayName(name);
//      expNames.add(displayName);
//      if (!displayName.equalsIgnoreCase(name)){
//        analysisResults.put(displayName, analysisResults.get(name));
//        analysisResults.remove(name);
//      }
//    }

    molBarCharts_.get(groupName).insertTab("Bar-chart", null, new BarChartPainter(BarChartPainter.TYPE_MOLECULE,groupName,moleculeName,analysisResults,analysisModule_.getExpNamesInSequence(),this,true,false,
        getMaxIsotopeForSetting(settingVO, maxIsotope,analysisResults), rtGrouped, false, settingVO, prefUnit, unit,
        analysisModule_.getCorrectionTypeISLookup().get(groupName), analysisModule_.getCorrectionTypeESLookup().get(groupName),new ArrayList<String>(analysisModule_.getModifications().get(groupName)),colorChooserDialog_),null,1);
    molBarCharts_.get(groupName).setSelectedIndex(1);
    return true; 
  }
  
  public boolean analyteGroupClicked(String moleculeName, String groupName, int maxIsotope, boolean rtGrouped, ResultDisplaySettingsVO settingVO, String prefUnit, String unit){
    molBarCharts_.get(groupName).remove(3);
    Hashtable<String,ResultCompVO> analysisResults = analysisModule_.getGroupedResults().get(groupName).get(moleculeName);
    Vector<String> groupNames = new Vector<String>(analysisModule_.getGroupNames());
    molBarCharts_.get(groupName).insertTab("Group bar-chart", null, new BarChartPainter(BarChartPainter.TYPE_MOLECULE,groupName,moleculeName,analysisResults,groupNames,this,false,true,
        getMaxIsotopeForSetting(settingVO, maxIsotope,analysisResults),rtGrouped, true, settingVO, prefUnit, unit,analysisModule_.getCorrectionTypeISLookup().get(groupName),
        analysisModule_.getCorrectionTypeESLookup().get(groupName),new ArrayList<String>(analysisModule_.getModifications().get(groupName)),colorChooserDialog_)
    ,null,3);
    molBarCharts_.get(groupName).setSelectedIndex(3);
    return true; 
  }
  
  public static int getMaxIsotopeForSetting(ResultDisplaySettingsVO settingVO, int maxIsotope,Hashtable<String,ResultCompVO> analysisResults){
    int maxAppliIsotope = maxIsotope+1;
    if (!settingVO.getType().equalsIgnoreCase(ResultDisplaySettingsVO.REL_MEASURED_CLASS_AMOUNT) && 
        !settingVO.getType().equalsIgnoreCase(ResultDisplaySettingsVO.REL_TOTAL_AMOUNT))
      maxAppliIsotope = StaticUtils.getMaxApplicableIsotope(analysisResults, maxIsotope);
    return maxAppliIsotope;
  }
  
  public boolean combinedAnalyteSelected(Vector<String> moleculeNames, String groupName, int maxIsotope, boolean rtGrouped, ResultDisplaySettingsVO settingVO, String prefUnit, String unit){
    molBarCharts_.get(groupName).remove(1);
    Hashtable<String,String> molNameHash = new Hashtable<String,String>();
    for (String molName : moleculeNames) molNameHash.put(molName, molName);
    Hashtable<String,Hashtable<String,ResultCompVO>> allAnalysisResults = new Hashtable<String,Hashtable<String,ResultCompVO>>(analysisModule_.getResults().get(groupName));
    Hashtable<String,Hashtable<String,ResultCompVO>> analysisResults = new Hashtable<String,Hashtable<String,ResultCompVO>>();
    Vector<String> expNames = new Vector<String>();
    Hashtable<String,String> expHash = new Hashtable<String,String>();
    for (String molName : allAnalysisResults.keySet()){
      if (molNameHash.containsKey(molName)){
        Hashtable<String,ResultCompVO> results = new Hashtable<String,ResultCompVO>();
        for (String name : analysisModule_.getExpNamesInSequence()){
//          String displayName  = getDisplayName(name);
//          if (!expHash.containsKey(displayName)){
            if (!expHash.containsKey(name)){
//            expNames.add(displayName);
//            expHash.put(displayName, displayName);
              expNames.add(name);
              expHash.put(name, name);
          }
//          results.put(displayName, allAnalysisResults.get(molName).get(name));
            results.put(name, allAnalysisResults.get(molName).get(name));
        }
        analysisResults.put(molName, results);
      }
    }
    molBarCharts_.get(groupName).insertTab("Bar-chart", null, new BarChartPainter(BarChartPainter.TYPE_MOLECULE,groupName,moleculeNames,analysisResults,expNames,this,true,false,
        StaticUtils.getMaxApplicableIsotopeHash(analysisResults, maxIsotope),rtGrouped,false,settingVO, prefUnit, unit, analysisModule_.getCorrectionTypeISLookup().get(groupName),
        analysisModule_.getCorrectionTypeESLookup().get(groupName),new ArrayList<String>(analysisModule_.getModifications().get(groupName)),colorChooserDialog_)
    ,null,1);
    molBarCharts_.get(groupName).setSelectedIndex(1);
    return true;
  }
  public boolean combinedAnalyteGroupSelected(Vector<String> moleculeNames, String groupName, int maxIsotope, boolean rtGrouped, ResultDisplaySettingsVO settingVO,
      String prefUnit, String unit){
    molBarCharts_.get(groupName).remove(3);
    Hashtable<String,String> molNameHash = new Hashtable<String,String>();
    Vector<String> groupNames = new Vector<String>(analysisModule_.getGroupNames());
    for (String molName : moleculeNames) molNameHash.put(molName, molName);
    Hashtable<String,Hashtable<String,ResultCompVO>> allAnalysisResults = new Hashtable<String,Hashtable<String,ResultCompVO>>(analysisModule_.getGroupedResults().get(groupName));
    Hashtable<String,Hashtable<String,ResultCompVO>> analysisResults = new Hashtable<String,Hashtable<String,ResultCompVO>>();
    for (String molName : allAnalysisResults.keySet()){
      if (molNameHash.containsKey(molName)){
        Hashtable<String,ResultCompVO> results = new Hashtable<String,ResultCompVO>();
        for (String name : groupNames){
          results.put(name, allAnalysisResults.get(molName).get(name));
        }
        analysisResults.put(molName, results);
      }
    }
    molBarCharts_.get(groupName).insertTab("Group bar-chart", null, new BarChartPainter(BarChartPainter.TYPE_MOLECULE,groupName,moleculeNames,analysisResults,groupNames,this,false,true,
        StaticUtils.getMaxApplicableIsotopeHash(analysisResults, maxIsotope),rtGrouped,true,settingVO, prefUnit, unit, analysisModule_.getCorrectionTypeISLookup().get(groupName),
        analysisModule_.getCorrectionTypeESLookup().get(groupName),new ArrayList<String>(analysisModule_.getModifications().get(groupName)),colorChooserDialog_)
    ,null,3);
    molBarCharts_.get(groupName).setSelectedIndex(3);
    return true;    
  }

  
  public boolean experimentClicked(String experimentName,String groupName, int maxIsotope, boolean rtGrouped, ResultDisplaySettingsVO settingVO, String prefUnit, String unit){
    molBarCharts_.get(groupName).remove(1);
    Hashtable<String,Hashtable<String,ResultCompVO>> analysisResults = analysisModule_.getResults().get(groupName);
    Hashtable<String,ResultCompVO> resultsForChart = new Hashtable<String,ResultCompVO>();
//    Vector<String> namesVector = analysisModule_.getAllMoleculeNames().get(groupName);
    Vector<String> namesVector = heatmaps_.get(groupName).getSelectedMoleculeNames();
    for (String moleculeName : analysisResults.keySet()){
      resultsForChart.put(moleculeName, analysisResults.get(moleculeName).get(experimentName));
    }
//    if (this.heatmaps_.get(groupName).hasNoInternalStandards()){
//      Hashtable<String,String> isNames = analysisModule_.getAllISNames().get(groupName);
//      for (String isName : isNames.keySet())
//        resultsForChart.remove(isName);
//      Vector<String> newNameVector = new Vector<String>();
//      for (String name: namesVector){
//        if (!isNames.containsKey(name))
//          newNameVector.add(name);
//      }
//      namesVector = new Vector<String>(newNameVector);
//    }
    molBarCharts_.get(groupName).insertTab("Bar chart", null,new BarChartPainter(BarChartPainter.TYPE_EXPERIMENT,groupName,getDisplayName(experimentName),resultsForChart,namesVector,this,false,false,
        getMaxIsotopeForSetting(settingVO, maxIsotope,resultsForChart),rtGrouped,false,settingVO, prefUnit, unit,analysisModule_.getCorrectionTypeISLookup().get(groupName),
        analysisModule_.getCorrectionTypeESLookup().get(groupName),new ArrayList<String>(analysisModule_.getModifications().get(groupName)),colorChooserDialog_)
    ,null,1);
    molBarCharts_.get(groupName).setSelectedIndex(1);
    return true;
  }
  
  public boolean experimentGroupClicked(String experimentGroupName,String groupName, int maxIsotope, boolean rtGrouped, ResultDisplaySettingsVO settingVO,
      String prefUnit, String unit){
    molBarCharts_.get(groupName).remove(3);
    Hashtable<String,Hashtable<String,ResultCompVO>> analysisResults = analysisModule_.getGroupedResults().get(groupName);
    Hashtable<String,ResultCompVO> resultsForChart = new Hashtable<String,ResultCompVO>();
//    Vector<String> namesVector = analysisModule_.getAllMoleculeNames().get(groupName);
    Vector<String> namesVector = groupHeatmaps_.get(groupName).getSelectedMoleculeNames();
    for (String moleculeName : analysisResults.keySet()){
      resultsForChart.put(moleculeName, analysisResults.get(moleculeName).get(experimentGroupName));
    }
//    if (this.groupHeatmaps_.get(groupName).hasNoInternalStandards()){
//      Hashtable<String,String> isNames = analysisModule_.getAllISNames().get(groupName);
//      for (String isName : isNames.keySet())
//        resultsForChart.remove(isName);
//      Vector<String> newNameVector = new Vector<String>();
//      for (String name: namesVector){
//        if (!isNames.containsKey(name))
//          newNameVector.add(name);
//      }
//      namesVector = new Vector<String>(newNameVector);
//    }
    molBarCharts_.get(groupName).insertTab("Group bar-chart", null,new BarChartPainter(BarChartPainter.TYPE_EXPERIMENT,groupName,getDisplayName(experimentGroupName),resultsForChart,namesVector,this,false,true,
        getMaxIsotopeForSetting(settingVO, maxIsotope,resultsForChart), rtGrouped,true,settingVO, prefUnit, unit, analysisModule_.getCorrectionTypeISLookup().get(groupName),
        analysisModule_.getCorrectionTypeESLookup().get(groupName),new ArrayList<String>(analysisModule_.getModifications().get(groupName)),colorChooserDialog_)
    ,null,3);
    molBarCharts_.get(groupName).setSelectedIndex(3);
    return true;
  }
  
  public String getDisplayName(String sampleName){
    if (groupDisplayNamesLookup_!=null && groupDisplayNamesLookup_.containsKey(sampleName)){
      return groupDisplayNamesLookup_.get(sampleName);
    }else{
      return this.expDisplayNamesLookup_.get(sampleName);
    }
  }
  
  public void setDisplayName(String sampleName, String displayName){
    this.expDisplayNamesLookup_.put(sampleName, displayName);
    for (String molGroup : heatmaps_.keySet()){
      heatmaps_.get(molGroup).generateHeatMap();
    }
    if (classOverviewPanel_!=null)
      classOverviewPanel_.refreshNames();
//    colorChooserDialog_.refreshNames();
  }
  
  public LinkedHashMap<String,String> getSampleResultFullPaths(){
    LinkedHashMap<String,String> fullPaths = new LinkedHashMap<String,String>();
    for (int i=0; i!=analysisModule_.getExpNamesInSequence().size(); i++) {
      String exp = analysisModule_.getExpNamesInSequence().get(i);
      fullPaths.put(exp, analysisModule_.getFullFilePath(exp).getAbsolutePath());
    }
    return fullPaths;
  }
  
  public void changeISStatus(String groupName, boolean isGrouped, boolean value){
    HeatMapDrawing heatMap = getCorrespodingHeatMap(groupName, isGrouped);
    if (heatMap !=null)
      heatMap.setISSelected(value);
  }
  
  private HeatMapDrawing getCorrespodingHeatMap(String groupName, boolean isGrouped){
    HeatMapDrawing heatMap = null;
    if (isGrouped)
      heatMap = heatmaps_.get(groupName);
    else{
      if (groupHeatmaps_.containsKey(groupName))
        heatMap = groupHeatmaps_.get(groupName);
    }
    return heatMap;
  }
  
  public void changeESStatus(String groupName, boolean isGrouped, boolean value){
    HeatMapDrawing heatMap = getCorrespodingHeatMap(groupName, isGrouped);
    if (heatMap !=null)
      heatMap.setESSelected(value);
  }
  
//  public void changeDoublePeakStatus(String groupName, boolean value){
//    HeatMapDrawing heatMap = getCorrespodingHeatMap(groupName, false);
//    if (heatMap !=null)
//      heatMap.setDoublePeakSelected(value);
//  }
  
  public void changeIsotopesUsed(String groupName, boolean isGrouped, int value){
    HeatMapDrawing heatMap = getCorrespodingHeatMap(groupName, isGrouped);
    if (heatMap !=null)
      heatMap.setSelectedIsotope(value);
  }
  
  @SuppressWarnings("unchecked")
  public void eliminateDoublePeaks(String groupName, String analyteName, String absFilePathStartEx, Vector<String> selectedMods, Vector<String> foundUpdateables){
    try{
      Hashtable<String,Boolean> showMods = new Hashtable<String,Boolean>();
      Hashtable<String,String> modHash = new Hashtable<String,String>();
      for (String modName : selectedMods) modHash.put(modName, modName);
      QuantificationResult result1 = LDAResultReader.readResultFile(absFilePathStartEx, showMods);
      Hashtable<String,LipidParameterSet> paramOfInterest = getParamByAnalyteName(analyteName, result1.getIdentifications().get(groupName), modHash);
      String[] nameAndRt = StaticUtils.extractMoleculeRtAndModFromMoleculeName(analyteName);
      if (nameAndRt[1]==null) nameAndRt[1] = "";
      String[] nameAndRtParam1 = StaticUtils.extractMoleculeRtAndModFromMoleculeName(analyteName);
      if (nameAndRtParam1[1]==null) nameAndRtParam1[1]="";
      if (paramOfInterest.size()>0){
        for (String updateablePath : foundUpdateables){
          Hashtable<String,Boolean> updateShowMods = new Hashtable<String,Boolean>();
          try{
            QuantificationResult result2 = LDAResultReader.readResultFile(updateablePath, updateShowMods);
            Vector<LipidParameterSet> updateParams = result2.getIdentifications().get(groupName);
            Hashtable<String,Vector<LipidParameterSet>> paramsToSelectOne = new Hashtable<String,Vector<LipidParameterSet>>();
            Vector<Integer> updateToRemove = new Vector<Integer>();
            Hashtable<String,Integer> positionToAdd = new Hashtable<String,Integer>();
            for (int j=0; j!=updateParams.size(); j++){
              LipidParameterSet param = updateParams.get(j);
              String[] nameAndRtParam2 = StaticUtils.extractMoleculeRtAndModFromMoleculeName(param.getNameString());
              if (nameAndRtParam2[1]==null) nameAndRtParam2[1]="";
              if (nameAndRtParam1[0].equalsIgnoreCase(nameAndRtParam2[0])){
                if (analysisModule_.neglectRtInformation(analyteName) || analysisModule_.isWithinRtGroupingBoundaries(Double.valueOf(nameAndRtParam1[1]), Double.valueOf(nameAndRtParam2[1]))){
                  if (modHash.containsKey(param.getModificationName())){
                    Vector<LipidParameterSet> paramsOfMod = new  Vector<LipidParameterSet>();
                    if (paramsToSelectOne.containsKey(param.getModificationName())) paramsOfMod = paramsToSelectOne.get(param.getModificationName());
                    paramsOfMod.add(param);
                    paramsToSelectOne.put(param.getModificationName(), paramsOfMod);
                    updateToRemove.add(0,j);
                    if (!positionToAdd.containsKey(param.getModificationName()))positionToAdd.put(param.getModificationName(), j);
                  }
                }
              }
            }
            for (int j : updateToRemove) updateParams.remove(j);
            Hashtable<String,LipidParameterSet> cleanedParams = new Hashtable<String,LipidParameterSet>();
            for (String modification : paramsToSelectOne.keySet()){
              Vector<LipidParameterSet> paramsOfMod = paramsToSelectOne.get(modification);
              LipidParameterSet paramOfInt = paramOfInterest.get(modification);
              LipidParameterSet remainingParam = null;
              for (LipidParameterSet set : paramsOfMod){
                if (remainingParam==null) remainingParam = set;
                else{
                  Vector<CgProbe> zeroProbes1 = remainingParam.getIsotopicProbes().get(0);
                  Vector<CgProbe> zeroProbes2 = set.getIsotopicProbes().get(0);
                  float refTime = paramOfInt.getIsotopicProbes().get(0).get(0).Peak;
                  float timeDifferenceMin = Float.MAX_VALUE;
                  for (CgProbe probe1 : zeroProbes1){
                    float timeDiff = StaticUtils.calculateDiff(refTime,probe1.Peak);
                    if (timeDiff<timeDifferenceMin)timeDifferenceMin = timeDiff;
                  }
                  boolean secondIsCloser = false;
                  for (CgProbe probe2 : zeroProbes2){
                    float timeDiff = StaticUtils.calculateDiff(refTime,probe2.Peak);
                    if (timeDiff<timeDifferenceMin){
                      timeDifferenceMin = timeDiff;
                      secondIsCloser = true;
                    }
                  }
                  if (secondIsCloser) remainingParam = set;
                }
              }
              Vector<Vector<CgProbe>> newIsotopicProbes = new Vector<Vector<CgProbe>>();
              for (int i=0; i!=remainingParam.getIsotopicProbes().size();i++){
                if (paramOfInt.getIsotopicProbes().size()>i){
                  if (paramOfInt.getIsotopicProbes().get(i).size()==1 && remainingParam.getIsotopicProbes().get(i).size()>1){
                    float retentionTime = paramOfInt.getIsotopicProbes().get(i).get(0).Peak;
                    float timeDifferenceMin = Float.MAX_VALUE;
                    CgProbe nearestProbe = null;
                    for (CgProbe aProbe : remainingParam.getIsotopicProbes().get(i)){
                      float timeDiff = aProbe.Peak-retentionTime;
                      if (timeDiff<0)
                        timeDiff = timeDiff*-1;
                      if (timeDiff<timeDifferenceMin){
                        nearestProbe = aProbe;
                        timeDifferenceMin = timeDiff;
                      }
                    }
                    Vector<CgProbe> newProbes = new Vector<CgProbe>();
                    newProbes.add(nearestProbe);
                    newIsotopicProbes.add(newProbes);
                  }else{
                    newIsotopicProbes.add(remainingParam.getIsotopicProbes().get(i));
                  }
                }
                if (i==0)remainingParam.setPreciseRt(newIsotopicProbes.get(0).get(0).Peak/60d);
              }
              float totalArea = 0f;
              Vector<CgProbe> probesVect = new Vector<CgProbe>();
              for (Vector<CgProbe> probes : newIsotopicProbes){
                for (CgProbe probe : probes){
                  totalArea+=probe.Area;
                  probesVect.add(probe);
                }  
              }
              remainingParam.setProbes(probesVect);
              remainingParam.Area = totalArea;
              cleanedParams.put(modification, remainingParam);
            }
            Vector<IntegerStringVO> positions = new Vector<IntegerStringVO>();
            for (String mod : positionToAdd.keySet())positions.add(new IntegerStringVO(mod,positionToAdd.get(mod)));
            Collections.sort(positions,new GeneralComparator("at.tugraz.genome.lda.vos.IntegerStringVO", "getValue", "java.lang.Integer"));
            int removalReduction = 0;
            for (IntegerStringVO position : positions){
              updateParams.add(position.getValue()-removalReduction,cleanedParams.get(position.getKey()));
              removalReduction = removalReduction+paramsToSelectOne.get(position.getKey()).size()-1;
            }
            try {
              QuantificationResultExporter.writeResultsToExcel(updateablePath, result2);
            }catch (Exception e) {e.printStackTrace();
            }
          } catch (ExcelInputFileException eif){
            //Comment: the graphical Warning message is shown in the readResultFile itself
          }
        }
        cleanupResultView();
        try
        {
        	acceptResultFiles();
        }
        catch (ExcelInputFileException ex)
        {
        	ex.printStackTrace();
        	new WarningMessage(new JFrame(), "Error", ex.getMessage());
        }
        for (int i=0; i!=resultTabs_.getTabCount();i++){
          if (resultTabs_.getTitleAt(i).equalsIgnoreCase(groupName))
            resultTabs_.setSelectedIndex(i);
        }      
      }else{
        new WarningMessage(new JFrame(), "Error", "There is something wrong with the file "+absFilePathStartEx+"! The "+analyteName+" in the group "+groupName+" does not exist!");
      }
    } catch (ExcelInputFileException eif){
      //Comment: the graphical Warning message is shown in the readResultFile itself
    }
  }
  
  public void addAnalytesEverywhereAtPosition(String groupName, Vector<String> analyteNames, Vector<String> absFilePathStartExps, Vector<String> selectedMods,
      Hashtable<String,Vector<AutoAnalyteAddVO>> updateablesAndAnalytesBefore, Vector<Integer> maxIsotopes, boolean exactProbePosition){
    //first, create a list of experimentPaths, so that the result files are not read more often than necessary
    //first key is the experiment, second key is the max MS-levels
    Hashtable<String,Integer> uniqueExampleExps = new Hashtable<String,Integer>();
    for (String exp : absFilePathStartExps) uniqueExampleExps.put(exp, 0);
    //second, create a list of experimentPaths of files that need to be updated
    Hashtable<String,String> uniqueUpdatePaths = new Hashtable<String,String>();
    for (Vector<AutoAnalyteAddVO> autoAnals : updateablesAndAnalytesBefore.values()){
      for (AutoAnalyteAddVO addVO : autoAnals) uniqueUpdatePaths.put(addVO.getResultFilePath(),addVO.getExperimentName());
    }
    try{
      //now read the template files and extract the corresponding template LipidParameterSets
      Hashtable<String,Boolean> showMods = new Hashtable<String,Boolean>();
      Hashtable<String,String> modHash = new Hashtable<String,String>();
      Hashtable<String,Hashtable<String,LipidParameterSet>> templateParams = new Hashtable<String,Hashtable<String,LipidParameterSet>>();
      for (String modName : selectedMods) modHash.put(modName, modName);
      for (String absFilePathStartEx : uniqueExampleExps.keySet()){
        //System.out.println("Now I am reading: "+absFilePathStartEx);
        QuantificationResult result1 = LDAResultReader.readResultFile(absFilePathStartEx, showMods);
        uniqueExampleExps.put(absFilePathStartEx, result1.getMsLevels().get(groupName));
        for (int i=0; i!=analyteNames.size(); i++){
          String analyteName = analyteNames.get(i);
          if (!absFilePathStartExps.get(i).equalsIgnoreCase(absFilePathStartEx))
            continue;
          Hashtable<String,LipidParameterSet> paramsOfInterest = getParamByAnalyteName(analyteName, result1.getIdentifications().get(groupName), modHash);
          //System.out.println("For "+analyteName+" I found "+paramsOfInterest.size()+" modifications.");
          templateParams.put(analyteName, paramsOfInterest);
        }
      }
      
      //now read, try to quantify, and update the corresponding experiments
      for (String updateAbsPath : uniqueUpdatePaths.keySet()){
        try{
          //System.out.println("file to be updated: "+updateAbsPath);
          QuantificationResult result2 = LDAResultReader.readResultFile(updateAbsPath, new Hashtable<String,Boolean>());
          String chromSetBasePath = StaticUtils.extractChromBaseName(updateAbsPath, uniqueUpdatePaths.get(updateAbsPath));
          String[] chromPaths = StringUtils.getChromFilePaths(chromSetBasePath+".chrom");
          LipidomicsAnalyzer analyzer = new LipidomicsAnalyzer(chromPaths[1],chromPaths[2],chromPaths[3],chromPaths[0],false);
          QuantificationThread.setAnalyzerProperties(analyzer);
          boolean hasRtInfo = QuantificationThread.hasRtInfo(result2.getIdentifications());
          boolean updateNecessary = false;
          //I have to start from the end of the analyte list, otherwise the position of the AddAnalyteAddVO is not correct any longer
          for (int i=(analyteNames.size()-1); i!=-1; i--){
            String analyteName = analyteNames.get(i);
            Hashtable<String,LipidParameterSet> paramsOfInterest = templateParams.get(analyteName);
            if (paramsOfInterest.size()==0)
              continue;
            //looking if there is an add
            AutoAnalyteAddVO addVO = null;
            for (AutoAnalyteAddVO voToAddBefore : updateablesAndAnalytesBefore.get(analyteName)){
              if (voToAddBefore.getResultFilePath().equalsIgnoreCase(updateAbsPath)){
                addVO = voToAddBefore;
                break;
              }
            }
            if (addVO==null)
              continue;
            //System.out.println("I am looking at: "+analyteName+" ; "+addVO.getPreviousElement());
            Vector<LipidParameterSet> updateParams = result2.getIdentifications().get(groupName);
            Vector<LipidParameterSet> updateParamsWOZeroAnalyte = new Vector<LipidParameterSet>();
            for (LipidParameterSet set: updateParams ){
              if (!set.getNameString().equalsIgnoreCase(analyteName) || set.Area>0f)
                updateParamsWOZeroAnalyte.add(set);
            }
            updateParams = updateParamsWOZeroAnalyte;
        
            int positionToAdd = 0;
            if (addVO.getPreviousElement()!=null && addVO.getPreviousElement().length()>0){
              String[] analyteBeforeNameAndRt = StaticUtils.extractMoleculeRtAndModFromMoleculeName(addVO.getPreviousElement());
              if (analyteBeforeNameAndRt[1]==null) analyteBeforeNameAndRt[1]="";
              for (int j=0;j!=updateParams.size();j++){
                String[] currentNameAndRt = StaticUtils.extractMoleculeRtAndModFromMoleculeName(updateParams.get(j).getNameString());
                if (currentNameAndRt[0].equalsIgnoreCase(analyteBeforeNameAndRt[0])){
                  if (analysisModule_.neglectRtInformation(analyteBeforeNameAndRt[0]) || analysisModule_.isWithinRtGroupingBoundaries(Double.valueOf(currentNameAndRt[1]), Double.valueOf(analyteBeforeNameAndRt[1]))){
                    positionToAdd = j+1;
//                    break;
                  }
                }
              }
            }
            //the position to add is known - now quantify each modification
            boolean isEmptyThere = true;
            for (String modName : paramsOfInterest.keySet()){
              LipidParameterSet templateParam = paramsOfInterest.get(modName);
              int charge = 1;
              if (templateParam.ProbeCount()>0)
                charge = templateParam.Probe(0).Charge;
              boolean oneThere = checkIsOneThere(analyteName,templateParam,updateParams);
//              if (oneThere)
//                System.out.println("There exists an empty one at :"+addVO.getExperimentName()+" ; "+modName);
              if (!oneThere){
                Vector<Vector<CgProbe>> probes = null;
                if (exactProbePosition)
                  probes = analyzer.calculatePeakAtExactProbePosition(templateParam,maxIsotopes.get(i),charge,uniqueExampleExps.get(absFilePathStartExps.get(i)));
                else
                  probes = analyzer.calculatePeakAtExactTimePosition(templateParam,maxIsotopes.get(i),charge,uniqueExampleExps.get(absFilePathStartExps.get(i)));
                // to calculate the total area
                  float totalArea = 0;
                  double rt = 0.0;
                  for (int k=0;k!=probes.size();k++){
                    Vector<CgProbe> isoProbes = probes.get(k);
                    Vector<Double> rts = new Vector<Double>();
                    for (CgProbe probe : isoProbes){
                      totalArea+=probe.Area;
                      rts.add((double)probe.Peak);
                    }
                    if (k==0 && hasRtInfo) rt = Calculator.mean(rts)/60d;
                  }
                  if (probes.size()==0 && hasRtInfo) rt = templateParam.getPreciseRT();
                  // INFO: Settings.emptyEntriesForQuantAnalNotFound() is responsible for creating empty entries if the analyte cannot be quantified                  
                  if (totalArea>0 || Settings.emptyEntriesForQuantAnalNotFound()){
                    LipidParameterSet paramToQuantify = new LipidParameterSet(templateParam.Mz[0], templateParam.getName(),
                      templateParam.getDoubleBonds(), templateParam.getModificationName(), rt, templateParam.getAnalyteFormula(),
                      templateParam.getModificationFormula(),templateParam.getCharge(), templateParam.getOhNumber());
                    
                    paramToQuantify.LowerMzBand = templateParam.LowerMzBand;
                    paramToQuantify.UpperMzBand = templateParam.UpperMzBand;
                    paramToQuantify.Area = totalArea;
                    paramToQuantify.setIsotopicProbes(probes);
                    String[] paramNameAndRt = StaticUtils.extractMoleculeRtAndModFromMoleculeName(paramToQuantify.getNameString());
                    boolean existsSameOne = false;
                    if (rt>0.0){
                      // check if the identified peak is already there
                      if ((positionToAdd-1)>-1){
                        String[] currentNameAndRt = StaticUtils.extractMoleculeRtAndModFromMoleculeName(updateParams.get(positionToAdd-1).getNameString());
                        if (currentNameAndRt[0].equalsIgnoreCase(paramNameAndRt[0]) && currentNameAndRt[1].equalsIgnoreCase(paramNameAndRt[1]))
                          existsSameOne = true;
                      }
                      for (int j=positionToAdd;j!=updateParams.size();j++){
                        String[] currentNameAndRt = StaticUtils.extractMoleculeRtAndModFromMoleculeName(updateParams.get(j).getNameString());
                        if (currentNameAndRt[0].equalsIgnoreCase(paramNameAndRt[0])){
                          if (paramToQuantify.getModificationName().equalsIgnoreCase(updateParams.get(j).getModificationName())){
                            if (currentNameAndRt[1].equalsIgnoreCase(paramNameAndRt[1])){
                              existsSameOne = true;
                              break;
                            }
                            else if (Double.valueOf(currentNameAndRt[1])<Double.valueOf(paramNameAndRt[1]))positionToAdd = j+1;
                            else{
                              positionToAdd = j;
                              break;
                            }  
                          }
                        }else break;                        
                      }
                    }
                    updateParams.add(positionToAdd, paramToQuantify);
                    if (rt>0.0)
                      positionToAdd++;
                    if (!existsSameOne)
                      isEmptyThere = false;
                 }
              }
            }
            if (!isEmptyThere){
              result2.getIdentifications().put(groupName, updateParams);
              //System.out.println("I found an update: "+analyteName+" ; "+addVO.getPreviousElement());
              updateNecessary = true;
            }
          }
          if (updateNecessary)
            QuantificationResultExporter.writeResultsToExcel(updateAbsPath, result2);
        } catch (Exception e) {
          // TODO Auto-generated catch block
          e.printStackTrace();
          new WarningMessage(new JFrame(), "Error", e.getMessage());
        }
      }
      cleanupResultView();
      try
      {
      	acceptResultFiles();
      }
      catch (ExcelInputFileException ex)
      {
      	ex.printStackTrace();
      	new WarningMessage(new JFrame(), "Error", ex.getMessage());
      }
      for (int i=0; i!=resultTabs_.getTabCount();i++){
        if (resultTabs_.getTitleAt(i).equalsIgnoreCase(groupName))
          resultTabs_.setSelectedIndex(i);
      }
    } catch (ExcelInputFileException eif){
      //Comment: the graphical Warning message is shown in the readResultFile itself
    }
  }
  

  public void eliminateAnalyteEverywhere(String groupName, Set<ResultCompVO> toRemove, Vector<String> selectedMods, Set<String> filePaths){
  	for (ResultCompVO compVO : toRemove)
  	{
  		String filePath = compVO.getAbsoluteFilePath();
  		filePaths.add(filePath);
  		ResultFileVO vo = analysisModule_.getResultFileVO(filePath);
  		if (vo == null) continue; //TODO: this should not happen, consider an error message.
  		
  		Vector<LipidParameterSet> currentParams = vo.getQuantificationResult().getIdentifications().get(groupName);
  		for (String mod : selectedMods)
  		{
  			Set<LipidParameterSet> params = compVO.getResultMolecule().getLipidParameterSets(mod);
  			currentParams.removeAll(params);
  		}
      vo.getQuantificationResult().getIdentifications().put(groupName, currentParams);
  	}
  	for (String filePath : filePaths)
  	{
  		try 
      {
        QuantificationResultExporter.writeResultsToExcel(filePath, analysisModule_.getResultFileVO(filePath).getQuantificationResult());
      }
      catch (Exception e) 
      {
        e.printStackTrace();
      }
  	}
    //TODO: here we could easily skip reading the files yet again and instead just update the heatmap of the one changed 'groupName'
    cleanupResultView();
    try
    {
    	acceptResultFiles();
    }
    catch (ExcelInputFileException ex)
    {
    	ex.printStackTrace();
    	new WarningMessage(new JFrame(), "Error", ex.getMessage());
    }
    for (int i=0; i!=resultTabs_.getTabCount();i++){
      if (resultTabs_.getTitleAt(i).equalsIgnoreCase(groupName))
        resultTabs_.setSelectedIndex(i);
    }    
  }
  
  public void addAnalyte(int position, AddAnalyteVO analyteDescrVO){
    QuantificationResult originalResult = new QuantificationResult(result_);
    this.currentSelected_ = position;
    Float exactMass = new Float(analyteDescrVO.getExactMass());
    LipidParameterSet set = new LipidParameterSet(exactMass, analyteDescrVO.getName(), 
        analyteDescrVO.getDoubleBonds(), analyteDescrVO.getModName(),
        Double.parseDouble(analyteDescrVO.getRt()), analyteDescrVO.getFormula(), analyteDescrVO.getModFormula(),
        new Integer(analyteDescrVO.getCharge()), analyteDescrVO.getOh());
      set.LowerMzBand = LipidomicsConstants.getCoarseChromMzTolerance(exactMass);
      set.UpperMzBand = LipidomicsConstants.getCoarseChromMzTolerance(exactMass);
      set.Area = 0;
      result_.getIdentifications().get(this.currentSelectedSheet_).add(position, set);
    try {
       storeResultsToExcel(position); 
    } catch (ExportException ex) {
      result_ = originalResult;
    }
  }
  
  public void removeAnalyte(int[] indices){
    QuantificationResult originalResult = new QuantificationResult(result_);
    List<Integer> inds = new ArrayList<Integer>();
    for (int ind:indices)inds.add(resultPositionToOriginalLoopkup_.get(ind));
    Collections.sort(inds);
    this.currentSelectedSheet_ =  (String)selectedSheet_.getSelectedItem();
    for (int i=(inds.size()-1); i>-1; i--){
      result_.getIdentifications().get(this.currentSelectedSheet_).remove(inds.get(i).intValue());
    }
    try {
      storeResultsToExcel(indices[0]);
    } catch (ExportException ex) {
      result_ = originalResult;
    }
  }
  
  public void removeMolecularSpecies(ArrayList<Integer> indices)
  {
  	Hashtable<Integer,LipidParameterSet> originalParams = new Hashtable<Integer,LipidParameterSet>();
  	for (int ind : indices)
  	{
  		LipidParameterSet param = getAnalyteInTableAtPosition(ind);
  		if (param instanceof LipidomicsMSnSet)
  		{
  			String molSpecies = ((LipidomicsTableModel)displayTable_.getModel()).getMSnIdentificationName(currentSelected_);
  			LipidomicsMSnSet paramMSn = (LipidomicsMSnSet)param;
  			try
  			{
  				LipidomicsMSnSet originalParam = new LipidomicsMSnSet(paramMSn);
  				originalParams.put(ind, originalParam); //storing a deep copy in case an ExportException occurs.
  				paramMSn.removeMolecularSpecies(molSpecies);
  			} 
  			catch (LipidCombinameEncodingException ex)
  			{
  				new WarningMessage(new JFrame(), "Error", 
  						String.format("The molecular species '%s' could not be deleted due to an internal error. Please contact the developers of the application if the problem persists.", molSpecies));
  				ex.printStackTrace();
  			}
  		}
  	}
  	if (!originalParams.isEmpty())
  	{
  		try
    	{
    		storeResultsToExcel(currentSelected_ > 0 ? (currentSelected_-1) : currentSelected_);
    	}
    	catch (ExportException ex) 
    	{
        for (int ind : indices) //resetting result_ to its original state.
        {
        	replaceResultParam(getAnalyteInTableAtPosition(ind), originalParams.get(ind));
        }
        new WarningMessage(new JFrame(), "Error", 
        		String.format("The following error occurred during the export: %s", ex.getMessage()));
      }
  	}
  }
  
  private void replaceResultParam(LipidParameterSet paramBefore, LipidParameterSet paramNew)
  {
  	int index = result_.getIdentifications().get(this.currentSelectedSheet_).indexOf(paramBefore);
  	result_.getIdentifications().get(currentSelectedSheet_).remove(index);
  	result_.getIdentifications().get(currentSelectedSheet_).add(index, paramNew);
  }

  public void updateColorChange()
  {
    for (String molName : molBarCharts_.keySet()){
      JTabbedPane tabbedPane = molBarCharts_.get(molName);
      for (int i=0;i!=tabbedPane.getTabCount();i++){
        if (tabbedPane.getComponentAt(i)!=null && tabbedPane.getComponentAt(i) instanceof BarChartPainter){
          BarChartPainter painter = (BarChartPainter)tabbedPane.getComponentAt(i);
          painter.repaint();
        }
      }
    }
    if (classOverviewPanel_!=null)
      classOverviewPanel_.repaint();
    if (classOverviewGroupPanel_!=null)
      classOverviewGroupPanel_.repaint();
  }
  
  public LipidParameterSet getAnalyteInTableAtPosition(int position){
    this.currentSelectedSheet_ = (String)selectedSheet_.getSelectedItem();
    return result_.getIdentifications().get(currentSelectedSheet_).get(resultPositionToOriginalLoopkup_.get(position));
  }
  
  public String getNameInTableAtPosition(int position)
  {
  	return ((LipidomicsTableModel)displayTable_.getModel()).getRowToName().get(currentSelected_);
  }
  
  private boolean checkIsOneThere(String analyteName, LipidParameterSet templateParam,Vector<LipidParameterSet> updateParams){
    String[] refNameAndRt =  StaticUtils.extractMoleculeRtAndModFromMoleculeName(analyteName);
    boolean oneThere = false;
    for (LipidParameterSet param : updateParams){
      String[] currentNameAndRt = StaticUtils.extractMoleculeRtAndModFromMoleculeName(param.getNameString());
      if (currentNameAndRt[0].equalsIgnoreCase(refNameAndRt[0]) && param.getModificationName().equalsIgnoreCase(templateParam.getModificationName())){
        if (analysisModule_.neglectRtInformation(analyteName) || analysisModule_.isWithinRtGroupingBoundaries(Double.valueOf(currentNameAndRt[1]), Double.valueOf(refNameAndRt[1]))){
          oneThere = true;
          break;
        }  
      }
    }
    return oneThere;
  }
  
  private void checkLicense()
  {
  	new LicenseChecker().checkLicense();
  }
  
  
  private class LicenseChangeListener implements ChangeListener{

    int currentSelectedIndex_;
    int lastSelectedIndex_;
   
    public LicenseChangeListener()
    {
      currentSelectedIndex_ = 0;
    }

    public void stateChanged(ChangeEvent e)
    {
    	if (LicenseChecker.isCheckLicense())
    	{
    		lastSelectedIndex_ = currentSelectedIndex_;
        currentSelectedIndex_ = mainTabs.getSelectedIndex();
        if (mainTabs.getSelectedIndex()==mainTabs.indexOfComponent(licensePanel_))
        {
        	LicenseChecker.showLicenseDialog();
          mainTabs.setSelectedIndex(lastSelectedIndex_);
        }
    	}
    }
  }
  
  private static Vector<File> sortFilesByName(Vector<File> filesToSort){
    Vector<File> newOrder = new Vector<File>(filesToSort.size());
    Hashtable<String,File> hash = new Hashtable<String,File>();
    for (File file : filesToSort) hash.put(file.getAbsolutePath(), file);
    Vector<String> pathVect = new Vector<String>(hash.keySet());
    Collections.sort(pathVect);
    for (String path : pathVect) newOrder.add(hash.get(path));
    return newOrder;
  }
  
    
 
  public void showMs2(int position){
    currentSelected_ = position;
    showMs2(null,false);
  }

  public boolean showMs2(LipidParameterSet set, boolean refreshL2dPainter){
    displayTolerancePanel_.getLockMzRange().setSelected(false);
    if (refreshL2dPainter)
      ms1ProbesWhileMs2Display_ = getAllProbesFromParams(set,null);
    else
      ms1ProbesWhileMs2Display_ = l2DPainter_.getAllSelectedProbes();      
    int charge = 1;
    
    displaysMs2_ = true;
    LipidParameterSet params = getAnalyteInTableAtPosition(currentSelected_);
    if (set !=null) params = set;
    if (params.ProbeCount()>0)
      charge = params.Probe(0).Charge;
    else if (params.getCharge()!=null&&params.getCharge()>1)
      charge = params.getCharge();
    float currentIsotopicMass = params.Mz[0];

    try {
      Hashtable<Integer,Vector<RangeColor>> rangeColors = StaticUtils.createRangeColorVOs(params,((LipidomicsTableModel)displayTable_.getModel()).getMSnIdentificationName(currentSelected_),
          result_.getFaHydroxyEncoding(), result_.getLcbHydroxyEncoding(), areTheseAlex123MsnFragments());
      int threeDMsLevel = 2;
      float tol = LipidomicsConstants.getMs2PrecursorTolerance(params.Mz[0]);
      Hashtable<Integer,Vector<String>> spectraRaw = reader_.getMsMsSpectra(params.Mz[0]-tol, params.Mz[0]+tol,-1f,-1f);

      float peakRt = 0f;
      if (params.getRt()!=null && params.getRt().length()>0)
        peakRt = Float.parseFloat(params.getRt())*60f;
      String[] chroms = null;
      int[] borders = null;
      int[] bordersThreeD = null;
      for (Integer level : spectraRaw.keySet()){
        String[] chromsLevel = reader_.translateSpectraToChroms(spectraRaw.get(level), level, LipidomicsConstants.getMs2ChromMultiplicationFactorForInt());
        int[] bordersLevel = reader_.getBordersOfLastMs2Extraction();
        if (level == threeDMsLevel){
          chroms = chromsLevel;
          bordersThreeD = bordersLevel;
        }
        if (borders==null){
          borders = bordersLevel;
        }else{
          if (bordersLevel[0]<borders[0]) borders[0] = bordersLevel[0];
          if (bordersLevel[2]<borders[2]) borders[2] = bordersLevel[2];
          if (bordersLevel[1]>borders[1]) borders[1] = bordersLevel[1];
          if (bordersLevel[3]>borders[3]) borders[3] = bordersLevel[3];
        }
      }
      if (chroms==null){
        throw new CgException("There are no MS2 spectra");
      }
      Hashtable<Integer,Float> allRetTimes = reader_.getMsmsRetentionTimes();
      Hashtable<Integer,Float> retTimes = reader_.getRetentionTimes(threeDMsLevel);
      float[] timeStartStop = getMSn3DTimeStartStopValues(retTimes,bordersThreeD[2],bordersThreeD[3],params.getIsotopicProbes().get(0));
      MSMapViewer viewer = MSMapViewerFactory.getMSMapViewer(chroms, retTimes,
          bordersThreeD[0],bordersThreeD[1],timeStartStop[0],timeStartStop[1],LipidomicsConstants.getMs2ChromMultiplicationFactorForInt(),1f,this,MSMapViewer.DISPLAY_TIME_MINUTES,true);
      if (rangeColors!=null)
        viewer.setMulticolorProbes(StaticUtils.get3DColorHash(rangeColors.get(threeDMsLevel)));
      viewer.setViewerSettings(true, true, true,LipidomicsConstants.getThreeDViewerMs2DefaultMZResolution(),LipidomicsConstants.getThreeDViewerMs2DefaultTimeResolution_());
      try{
        viewer.init();
      }catch (Exception cgx) {
        cgx.printStackTrace();
        @SuppressWarnings("unused")
        WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "The 3D Viewer cannot be started: "+cgx.getMessage());
      }
      viewer.removeSaveLipidomicsSettings();
      @SuppressWarnings("rawtypes")
      Vector<Hashtable> rtNrSpectraAndPrecursor = reader_.getRtNrSpectrumHash(spectraRaw);
      @SuppressWarnings("unchecked")
      Hashtable<Integer,String> scanNrSpectrumHash = (Hashtable<Integer,String>)rtNrSpectraAndPrecursor.get(0);
      @SuppressWarnings("unchecked")
      Hashtable<Integer,Vector<Double>> scanNrPrecursorHash = (Hashtable<Integer,Vector<Double>>)rtNrSpectraAndPrecursor.get(1);
      @SuppressWarnings("unchecked")
      Hashtable<Integer,Integer> scanNrLevelHash = (Hashtable<Integer,Integer>)rtNrSpectraAndPrecursor.get(2);
      if (this.displayTolerancePanel_.getShow2D().isSelected()){
        int extendedStart = (borders[0]*9)/10;
        int extendedStop = (borders[1]*21)/20;
        Lipidomics2DSpectraChromPainter spectrumPainter = new Lipidomics2DSpectraChromPainter(analyzer_,scanNrSpectrumHash, scanNrPrecursorHash, scanNrLevelHash, allRetTimes,
            peakRt, extendedStart,extendedStop,LipidomicsConstants.getMs2ChromMultiplicationFactorForInt(),tol*2,this,
            extendedStart,extendedStop, true,new Vector<CgProbe>(),new Vector<CgProbe>(),1,1, this.relAbund_.isSelected(),
            params,rangeColors, new Double(annotationThreshold_.getText()),shotgunIsDisplayed_);

        this.spectrumSelected_.setText(spectrumPainter.getSpectSelectedText());
        this.rtSelected_.setText(spectrumPainter.getRtSelectedText());
        displayPrecursorMasses(spectrumPainter.getPrecursorMassSelected());
        this.msLevelSelected_.setText(spectrumPainter.getMsLevelSelected());
        float[] range = spectrumPainter.getRTRange(); 
        viewer.setCurrent2DTimeRange(range[0], range[1]);       
        makeDisplayRemoveOperations(false);
        this.viewer_ = viewer;

        
        majorSplitPane_ = new JSplitPane(JSplitPane.VERTICAL_SPLIT);
        majorSplitPane_.setDividerSize(2);
        JPanel topSplitPanel = new JPanel();
        topSplitPanel.setLayout(new BorderLayout());
        topSplitPane_ = new JSplitPane(JSplitPane.VERTICAL_SPLIT);
        topSplitPanel.add(topSplitPane_,BorderLayout.CENTER);
        majorSplitPane_.setTopComponent(topSplitPanel);
        topSplitPane_.setDividerSize(2);
        topSplitPane_.setTopComponent(viewer_);
        topSplitPane_.setBottomComponent(this.l2dPanel_);
        //topSplitPane_.setDividerLocation(0.6d);
        
        spectrumPainter_ = spectrumPainter;
        spectrumPanel_.add(spectrumPainter_,BorderLayout.CENTER);
        majorSplitPane_.setBottomComponent(this.spectrumPanel_);
        spectrumPainter.setBackground(Color.WHITE);
      
        this.displayPanel_.add(majorSplitPane_,BorderLayout.CENTER);
        this.displayPanel_.validate();
        majorSplitPane_.setDividerLocation(0.6d);
      
        
        spectrumPainter_.draw2DDiagram(extendedStart,extendedStop, true);
        
        
        // this is for setting the displayed chromatogram to the isotope 0
        if (refreshL2dPainter || isotope_.getSelectedIndex()!=0){
          isotope_.setSelectedIndex(0);
          float startFloat = currentIsotopicMass-Float.parseFloat(this.displayTolerancePanel_.getDisplayMinusTolerance().getText());
          float stopFloat = currentIsotopicMass+Float.parseFloat(this.displayTolerancePanel_.getDisplayPlusTolerance().getText());
          Vector<CgProbe> storedProbes = ms1ProbesWhileMs2Display_.get(0);
          Vector<CgProbe> selectedProbes = ms1ProbesWhileMs2Display_.get(1);
          if (l2DPainter_!=null){
            l2DPainter_.getGraphics().dispose();
            l2dPanel_.remove(l2DPainter_);
            l2DPainter_ = null;
          }
          
          try {
            String[] rawLines = reader_.getRawLines(startFloat, stopFloat, result_.getMsLevels().get(currentSelectedSheet_));
            Hashtable<Integer,Float> rtTimes = reader_.getRetentionTimes();
            Lipidomics2DPainter l2DPainter = new Lipidomics2DPainter(analyzer_,rawLines, reader_.getRetentionTimesOriginal(), rtTimes,
                startFloat,stopFloat,reader_.getMultiplicationFactorForInt_()/reader_.getLowestResolution_(),params_.LowerMzBand*2,this,
                MSMapViewer.DISPLAY_TIME_MINUTES,currentIsotopicMass-params.LowerMzBand, currentIsotopicMass+params.UpperMzBand, false,
                storedProbes,selectedProbes,Integer.parseInt((String)this.isotope_.getSelectedItem()),charge,result_.getMsLevels().get(currentSelectedSheet_),
                shotgunIsDisplayed_);
            
            l2DPainter.preChromatogramExtraxtion(currentIsotopicMass-params.LowerMzBand, currentIsotopicMass+params.UpperMzBand);
            l2DPainter_ = l2DPainter;
            l2DPainter_.setBackground(Color.WHITE);
            l2dPanel_.add(l2DPainter,BorderLayout.CENTER);
          }
          catch (CgException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
          }

        }
        l2DPainter_.repaint();
        l2DPainter_.paintMs2Position(Float.parseFloat(rtSelected_.getText()));
        majorSplitPane_.repaint();
        topSplitPane_.setResizeWeight(0.58);
        showMs2InputElements();
      }else{
        this.makeDisplayRemoveOperations();
        this.viewer_ = viewer;
        this.displayPanel_.add(viewer_,BorderLayout.CENTER);
        this.displayPanel_.validate();
      }
      

      
    }catch (CgException cgx) {
      new WarningMessage(new JFrame(), "Error", "The MS/MS cannot be displayed: "+cgx.getMessage());
      this.displaysMs2_ = false;
      return false;
    }catch (LipidCombinameEncodingException cgx) {
      new WarningMessage(new JFrame(), "Error", "The molecular species name cannot be decoded: "+cgx.getMessage());
      this.displaysMs2_ = false;
      return false;

    }catch (Exception cgx) {
      cgx.printStackTrace();
      new WarningMessage(new JFrame(), "Warning", "The 3D Viewer cannot be started: "+cgx.getMessage());
    }
    return true;
  }
  
  private Vector<Vector<CgProbe>> getAllProbesFromParams(LipidParameterSet params,Vector<Vector<CgProbe>> previouslyselectedProbes){
    Vector<Vector<CgProbe>> allProbes = new Vector<Vector<CgProbe>>(2);
    Vector<CgProbe> storedProbes = new Vector<CgProbe>();
    Vector<CgProbe> selectedProbes = new Vector<CgProbe>();
    for (int i=0; i!=params.ProbeCount();i++){
      storedProbes.add(params.Probe(i));
    }
    if (previouslyselectedProbes!=null){
      storedProbes = previouslyselectedProbes.get(0);
      selectedProbes = previouslyselectedProbes.get(1);
    }
    allProbes.add(storedProbes);
    allProbes.add(selectedProbes);
    return allProbes;
  }
  
  public void exportMzTab(File exportFile, short speciesTypeParam, boolean exportDoubleBondPositions)
  {
    
	    Metadata metadata = new Metadata();
    metadata.setIdConfidenceMeasure(new ArrayList<Parameter>());
    //TODO: this has to be generated somehow differently
	    metadata.setMzTabID("1");
	    Software software = new Software();
	    software.setId(1);
	    software.setParameter(new Parameter().name("LipidDataAnalyzer").value(Settings.VERSION));
	    metadata.addSoftwareItem(software);
	    Parameter quantMethod = new Parameter().cvLabel("MS").cvAccession("MS:1002019").name("label free raw feature quantitation").value(null);
    metadata.setQuantificationMethod(quantMethod);
	    List<CV> cvs = new ArrayList<CV>();
	    CV label = new CV().id(1).label("MS").fullName("PSI-MS controlled vocabulary").version("20-06-2018").uri("https://www.ebi.ac.uk/ols/ontologies/ms");
	    cvs.add(label);
    label = new CV().id(2).label("PRIDE").fullName("PRIDE PRoteomics IDEntifications (PRIDE) database controlled vocabulary").version("14-06-2018").uri("https://www.ebi.ac.uk/ols/ontologies/pride");
    cvs.add(label);
    
    int cvId = 3;
    if (LipidomicsConstants.isCvSepRequired()){
      label = new CV().id(cvId).label("SEP").fullName("Sample Processing and Separation Techniques Ontology").version("1.070708").uri("http://purl.bioontology.org/ontology/SEP");
      cvs.add(label);
      cvId++;
    }
    if (LipidomicsConstants.isCvChmoRequired()){
      label = new CV().id(cvId).label("CHMO").fullName("Chemical Methods Ontology").version("14-11-2019").uri("http://purl.obolibrary.org/obo/chmo.owl");
      cvs.add(label);
      cvId++;
    }
    if (LipidomicsConstants.isCvNcbiTaxonRequired()){      
      label = new CV().id(cvId).label("NCBITaxon").fullName("NCBI organismal classification").version("2018-03-02").uri("https://www.ebi.ac.uk/ols/ontologies/ncbitaxon");
      cvs.add(label);
      cvId++;
    }
    if (LipidomicsConstants.isClRequired()){      
      label = new CV().id(cvId).label("CL").fullName("Cell Ontology").version("2018-07-07").uri("https://www.ebi.ac.uk/ols/ontologies/cl");
      cvs.add(label);
      cvId++;
    }
    if (LipidomicsConstants.isCvBtoRequired()){      
      label = new CV().id(cvId).label("BTO").fullName("BRENDA tissue / enzyme source").version("2016-05-05").uri("https://www.ebi.ac.uk/ols/ontologies/bto");
      cvs.add(label);
      cvId++;
    }
    if (LipidomicsConstants.isCvDoidRequired()){      
      label = new CV().id(cvId).label("DOID").fullName("Human Disease Ontology").version("2018-07-05").uri("https://www.ebi.ac.uk/ols/ontologies/doid");
      cvs.add(label);
      cvId++;
    } 
    metadata.setCv(cvs);
    
	    List<Database> databases = new ArrayList<Database>();
	    //TODO: exchange with corresponding MS accession
	////    Database database = new Database().id(1).param(new Parameter().cvLabel("MS").cvAccession("MS:XXXXXXX").name("Lipid Data Analyzer software for lipid quantification"));
	    Database database = new Database().id(1).param(new Parameter().name("LipidDataAnalyzer2").value("lda2"));
	    database.setPrefix("lda2");
	    database.setVersion(Settings.VERSION);
	    database.setUri("https://github.com/ThallingerLab/LDA2");
	    databases.add(database);
	    metadata.setDatabase(databases);
	    
	    List<Instrument> instruments = new ArrayList<Instrument>();
	    Instrument instrument = LipidomicsConstants.getMzTabInstrument();
	    instruments.add(instrument);
	    if (instrument!=null)
	      metadata.setInstrument(instruments);

	    for (Contact contact : LipidomicsConstants.getMzTabContacts()){
	      metadata.addContactItem(contact);
	    }
	    List<SampleProcessing> processings = LipidomicsConstants.getMzTabSampleprocessings();
	    metadata.setSampleProcessing(processings);
	    metadata.setPublication(LipidomicsConstants.getMzTabPublications());
	    if (LipidomicsConstants.getMzTabSample()!=null){
	      List<Sample> samples = new ArrayList<Sample>();
	      samples.add(LipidomicsConstants.getMzTabSample());
	      metadata.setSample(samples);
	    }
	    metadata.setSmallMoleculeQuantificationUnit(new Parameter().cvLabel("PRIDE").cvAccession("PRIDE:0000330").name("Arbitrary quantification unit"));
	    metadata.setSmallMoleculeFeatureQuantificationUnit(new Parameter().cvLabel("PRIDE").cvAccession("PRIDE:0000330").name("Arbitrary quantification unit"));
	    //TODO: exchange with corresponding MS accession
	    metadata.setSmallMoleculeIdentificationReliability(new Parameter().cvLabel("MS").cvAccession("MS:1002896").name("compound identification confidence level"));
    List<Parameter> idConfidenceMeasures = new ArrayList<Parameter>();
    idConfidenceMeasures.add(new Parameter().id(1).cvLabel("MS").cvAccession("MS:1002890").name("fragmentation score"));
    metadata.setIdConfidenceMeasure(idConfidenceMeasures);
	    
    try {
	      Hashtable<String,Integer> expToMsRun = new Hashtable<String,Integer>();
	      resultsShowModification_ = new Hashtable<String,Boolean>();
    
	      Hashtable<String,QuantificationResult> originalExcelResults = new Hashtable<String,QuantificationResult>();
	      Hashtable<String,MsRun> msRuns = new Hashtable<String,MsRun>();
	      LinkedHashMap<String,String> fullExpPaths = getSampleResultFullPaths();
	      boolean isotopicDistributionChecked = true;
	      Hashtable<String,Short> polarities = new Hashtable<String,Short>();
	      HydroxyEncoding faEncoding = null;
	      HydroxyEncoding lcbEncoding = null;
	      for (int i=0; i!=analysisModule_.getExpNamesInSequence().size(); i++) {
	        String exp = analysisModule_.getExpNamesInSequence().get(i);
	        String chromFileBase = StaticUtils.extractChromBaseName(fullExpPaths.get(exp),exp);
	        MsRun run = new MsRun();
	        run.setId((i+1));
	        run.setLocation("file://"+chromFileBase.replaceAll("\\\\", "/").replaceAll(" ", "%20")+".chrom");       
	        run.setFragmentationMethod(LipidomicsConstants.getFragmentationMethods());
	        run.setFormat(new Parameter().cvLabel("MS").cvAccession("MS:1002966").name("The Lipid Data Analyzer native chrom format."));
	        run.setIdFormat(new Parameter().cvLabel("MS").cvAccession("MS:1000776").name("scan number only nativeID format"));
	        run.setHashMethod(new Parameter().cvLabel("MS").cvAccession("MS:1000569").name("SHA-1"));
	        
	        Assay assay = new Assay();
	        assay.setId(i+1);
	        assay.addMsRunRefItem(run);
	        msRuns.put(exp, run);
////      assay.setQuantificationReagent(quantificationMethod);
	        metadata.addAssayItem(assay);
	        expToMsRun.put(exp, (i+1));
	        metadata.addMsRunItem(run);
	        originalExcelResults.put(exp, LDAResultReader.readResultFile(fullExpPaths.get(exp), new Hashtable<String,Boolean>()));
	        if (originalExcelResults.get(exp).getFaHydroxyEncoding()!=null) faEncoding = originalExcelResults.get(exp).getFaHydroxyEncoding();
	        if (originalExcelResults.get(exp).getLcbHydroxyEncoding()!=null) lcbEncoding = originalExcelResults.get(exp).getLcbHydroxyEncoding();
	        if (!originalExcelResults.get(exp).getConstants().getRespectIsotopicDistribution())
	          isotopicDistributionChecked = false;
	        polarities.put(exp, SmallMztabMolecule.POLARITY_UNKNOWN);
////	      factory.addAbundanceOptionalColumn(assay);
	      }
	      if (isotopicDistributionChecked)
	        software.addSettingItem("isotope_pattern_checked");
	      
	      Hashtable<String,LinkedHashMap<String,Boolean>> adductsSortedByAbundance = LDAExporter.extractAdductsSortedByAbundance(this.heatmaps_.keySet(),originalExcelResults);
      //check if the settings type is everywhere the same, otherwise, use ResultDisplaySettingsVO.REL_VALUE
	      String valueType = null;
	      boolean allTheSame = true;
	      for (String molGroup : heatmaps_.keySet()){
	        String type = heatmaps_.get(molGroup).getValueType();
	        if (valueType==null)
	          valueType = type;
	        else if (!type.equalsIgnoreCase(valueType))
	          allTheSame = false;
      }
	      if (!allTheSame)
	        valueType = ResultDisplaySettingsVO.REL_VALUE;
	      LinkedHashMap<String,Vector<String>> expsOfGroup = getSamplesOfGroups();
	      int count = 0;
	      for (String group : expsOfGroup.keySet()){
        StudyVariable studyVariable = new StudyVariable();
        studyVariable.setId(count+1);
        studyVariable.setName(group);
        studyVariable.setDescription(group);
        studyVariable.setAverageFunction(new Parameter().cvLabel("MS").cvAccession("MS:1002962").name("The arithmetic mean"));
        studyVariable.setVariationFunction(new Parameter().cvLabel("MS").cvAccession("MS:1002963").name("The coefficient of variation"));
        for (String exp : expsOfGroup.get(group)){
          Integer assayId = expToMsRun.get(exp);
          for (Assay assay : metadata.getAssay()){
            if (assay.getId() == assayId)
              studyVariable.addAssayRefsItem(assay);
          }
        }
        metadata.addStudyVariableItem(studyVariable);
	        count++;
      }  
	  
	      MzTab mzTabFile = new MzTab();
	      mzTabFile.setMetadata(metadata);
      Hashtable<String,Vector<String>> correctAnalyteSequence = this.analysisModule_.getAllMoleculeNames();
      int summaryId = 1;
      int featureId = 1;
      int evidenceId = 1;
      int evidenceGroupingId = 1;
      short polarity;
      boolean containsSMESection = false;
	      for (String molGroup:this.heatmaps_.keySet()) {
	        boolean exportDoubleBondPositionsForClass = exportDoubleBondPositions; 
	        short speciesType = speciesTypeParam;
	        //when lipid species level is selected, export double bond positions only for lipid species consisting of one FA chain
	        if (exportDoubleBondPositions) {
	          if (analysisModule_.getNrOfChainsOfClass().get(molGroup) > 1 && speciesType == LipidomicsConstants.EXPORT_ANALYTE_TYPE_SPECIES) {
	            exportDoubleBondPositionsForClass = false;
	          } else if (analysisModule_.getNrOfChainsOfClass().get(molGroup) == 1 && speciesType != LipidomicsConstants.EXPORT_ANALYTE_TYPE_SPECIES) {
	            speciesType = LipidomicsConstants.EXPORT_ANALYTE_TYPE_SPECIES;
	          }
	        }
	        Hashtable<String,String> selectedMolHash = new Hashtable<String,String>();
	        HeatMapDrawing heatmap = this.heatmaps_.get(molGroup);
	//        HeatMapDrawing groupHeatmap = null;
	//        if (analysisModule_.getGroupNames()!=null && analysisModule_.getGroupNames().size()>0) groupHeatmap = this.groupHeatmaps_.get(molGroup);
	        for (String molName: heatmap.getSelectedMoleculeNames()){
	          selectedMolHash.put(molName, molName);
	        }
	        Hashtable<String,Hashtable<String,Vector<Double>>> results = heatmap.getResultValues(valueType);
	//        Hashtable<String,Hashtable<String,Vector<Double>>> resultsGroup = null;
	//        if (groupHeatmap!=null) resultsGroup = groupHeatmap.getResultValues(valueType);
        int maxIsotopes = heatmap.getSelectedIsotope();
        LinkedHashMap<String,Boolean> adductsSorted = adductsSortedByAbundance.get(molGroup);
        for (String molName: correctAnalyteSequence.get(molGroup)){
          if (selectedMolHash.containsKey(molName)){
            Hashtable<String,Vector<Double>> resultsMol = results.get(molName); 
//            Hashtable<String,Vector<Double>> resultsGroupMol = null;
//	            if (resultsGroup != null) resultsGroupMol = resultsGroup.get(molName);
//	            
          
	//            String identifier = null;
	//              identifier = "IS_" + molGroup + molName.substring(internalStandardPref.length());
	 //             identifier = "ES_" + molGroup + molName.substring(externalStandardPref.length());
	//            else
	//              identifier = molGroup + molName;
            SmallMztabMolecule molecule = MztabUtils.createSmallMztabMolecule(speciesType, exportDoubleBondPositionsForClass, 
                summaryId, featureId, evidenceId,evidenceGroupingId, maxIsotopes,analysisModule_, msRuns, originalExcelResults, 
                molGroup, molName, resultsMol,adductsSorted, expsOfGroup, faEncoding, lcbEncoding);
	            if (molecule==null)
	              continue;
	            summaryId = molecule.getCurrentSummaryId();
            featureId = molecule.getCurrentFeatureId();
            evidenceId = molecule.getCurrentEvidenceId();
            evidenceGroupingId = molecule.getCurrentEvGroupingId();
            for (SmallMoleculeSummary summary : molecule.getSummary()){
		              mzTabFile.addSmallMoleculeSummaryItem(summary);
	            }
	            for (SmallMoleculeFeature feature : molecule.getFeatures()){
	              mzTabFile.addSmallMoleculeFeatureItem(feature);
	            }
	            for (SmallMoleculeEvidence evidence : molecule.getEvidence()){
	              mzTabFile.addSmallMoleculeEvidenceItem(evidence);
	              if (!containsSMESection)
	                containsSMESection = true;
	            }
	            for (int i=0; i!=analysisModule_.getExpNamesInSequence().size(); i++) {
	              String exp = analysisModule_.getExpNamesInSequence().get(i);
	              polarity = molecule.getPolarity().get(exp);
	              if (polarity==SmallMztabMolecule.POLARITY_UNKNOWN || polarity==polarities.get(exp))
	                continue;
	              if (polarity==SmallMztabMolecule.POLARITY_BOTH)
	                polarities.put(exp, SmallMztabMolecule.POLARITY_BOTH);
	              else if (polarity==SmallMztabMolecule.POLARITY_POSITIVE)
	                polarities.put(exp, polarities.get(exp)==SmallMztabMolecule.POLARITY_NEGATIVE ? SmallMztabMolecule.POLARITY_BOTH : SmallMztabMolecule.POLARITY_POSITIVE);
	              else if (polarity==SmallMztabMolecule.POLARITY_NEGATIVE)
	                polarities.put(exp, polarities.get(exp)==SmallMztabMolecule.POLARITY_POSITIVE ? SmallMztabMolecule.POLARITY_BOTH : SmallMztabMolecule.POLARITY_NEGATIVE);
	            }
          }
	        }
	      }
	      
	      
	      for (String exp : analysisModule_.getExpNamesInSequence()) {
	        List<Parameter> polaries = new ArrayList<Parameter>();
	        if (polarities.get(exp) == SmallMztabMolecule.POLARITY_POSITIVE || polarities.get(exp) == SmallMztabMolecule.POLARITY_BOTH)
	          polaries.add(new Parameter().cvLabel("MS").cvAccession("MS:1000130").name("positive scan"));
        if (polarities.get(exp) == SmallMztabMolecule.POLARITY_NEGATIVE || polarities.get(exp) == SmallMztabMolecule.POLARITY_BOTH)
          polaries.add(new Parameter().cvLabel("MS").cvAccession("MS:1000129").name("negative scan"));
        msRuns.get(exp).setScanPolarity(polaries);
	      }
      
	      //BufferedOutputStream stream = new BufferedOutputStream(new FileOutputStream(exportFile));
	      MzTabWriter<List<ValidationMessage>> writer = null;
	      if (containsSMESection)
	        writer = new MzTabValidatingWriter();
	      else
	        writer = new MzTabValidatingWriter(new MzTabValidatingWriter.WriteAndParseValidator(System.out, Level.Warn, 100),new MzTabWriterDefaults(), true);

	      OutputStreamWriter outwriter = new OutputStreamWriter(new FileOutputStream(exportFile.getAbsolutePath()), StandardCharsets.UTF_8);
	      writer.write(outwriter, mzTabFile);
	      outwriter.close();
	    }
	    catch (IOException e) {
	      e.printStackTrace();
	    }
	    catch (NumberFormatException e) {
	      e.printStackTrace();
	    }
	    catch (CalculationNotPossibleException e) {
	      e.printStackTrace();
	    }
    catch (ExcelInputFileException e) {
	      e.printStackTrace();
	    }
    catch (ExportException | LipidCombinameEncodingException e ) {
      e.printStackTrace();
    }
    catch (SpectrummillParserException e) {
      e.printStackTrace();
      new WarningMessage(new JFrame(), "Error", e.getMessage());
    }
    catch (RetentionTimeGroupingException e) {
      new WarningMessage(new JFrame(), "Error", e.getMessage());
    }
	}
  
  public void exportRdb(File exportFile){
    System.out.println("RDB export");
    // save the internal and external standard prefix
    String internalStandardPref = internalStandardSelection_.getText();
    String externalStandardPref = externalStandardSelection_.getText();

    
    try {
      Hashtable<String,Hashtable<String,String>> acceptedMolecules = new Hashtable<String,Hashtable<String,String>>();
      LinkedHashMap<String,Integer> classSequence = new LinkedHashMap<String,Integer>();
      Hashtable<String,Integer> maxIsotopes = new Hashtable<String,Integer>();
      for (String molGroup:this.heatmaps_.keySet()) {
        classSequence.put(molGroup, 1);
        Hashtable<String,String> molsOfClass = new Hashtable<String,String>();
        HeatMapDrawing heatmap = this.heatmaps_.get(molGroup);
        for (String molName: heatmap.getSelectedMoleculeNames()){
          molsOfClass.put(molName, molName);
        }
        acceptedMolecules.put(molGroup, molsOfClass);
        maxIsotopes.put(molGroup, heatmap.getSelectedIsotope());
      }
      RdbOutputWriter rdbWriter = new RdbOutputWriter(internalStandardPref,externalStandardPref);
      rdbWriter.write(exportFile.getAbsolutePath(), analysisModule_, classSequence, null, acceptedMolecules,maxIsotopes,null,null,false);
    }
    catch (ExcelInputFileException e) {
      e.printStackTrace();
      @SuppressWarnings("unused")
      WarningMessage dlg = new WarningMessage(new JFrame(), "Error", "The Alex123-file cannot be written: "+e.getMessage());
    }
    catch (RdbWriterException rwe) {
      rwe.printStackTrace();
      @SuppressWarnings("unused")
      WarningMessage dlg = new WarningMessage(new JFrame(), "Error", "The Alex123-file cannot be written: "+rwe.getMessage());
    }
    catch (ChemicalFormulaException | LipidCombinameEncodingException cfe) {
      cfe.printStackTrace();
      @SuppressWarnings("unused")
      WarningMessage dlg = new WarningMessage(new JFrame(), "Error", "The Alex123-file cannot be written: "+cfe.getMessage());
    }
  }
  
  
  @SuppressWarnings("unchecked")
  public void exportMaf(File exportFile){
    try {
      BufferedOutputStream stream = new BufferedOutputStream(new FileOutputStream(exportFile));
      Hashtable<String,Hashtable<String,Hashtable<String,Vector<LipidomicsMSnSet>>>> results = new Hashtable<String,Hashtable<String,Hashtable<String,Vector<LipidomicsMSnSet>>>>();
      String headerLine = "";
      headerLine += "\"database_identifier\"";
      headerLine += "\t"+"\"chemical_formula\"";
      headerLine += "\t"+"\"smiles\"";
      headerLine += "\t"+"\"inchi\"";
      headerLine += "\t"+"\"metabolite_identification\"";
      headerLine += "\t"+"\"mass_to_charge\"";
      headerLine += "\t"+"\"fragmentation\"";
      headerLine += "\t"+"\"modifications\"";
      headerLine += "\t"+"\"charge\"";
      headerLine += "\t"+"\"retention_time\"";
      headerLine += "\t"+"\"taxid\"";
      headerLine += "\t"+"\"species\"";
      headerLine += "\t"+"\"database\"";
      headerLine += "\t"+"\"database_version\"";
      headerLine += "\t"+"\"reliability\"";
      headerLine += "\t"+"\"uri\"";
      headerLine += "\t"+"\"search_engine\"";
      headerLine += "\t"+"\"search_engine_score\"";
      headerLine += "\t"+"\"smallmolecule_abundance_sub\"";
      headerLine += "\t"+"\"smallmolecule_abundance_stdev_sub\"";
      headerLine += "\t"+"\"smallmolecule_abundance_std_error_sub\"";

      for (int i=0; i!=analysisModule_.getExpNamesInSequence().size(); i++) {
        String exp = analysisModule_.getExpNamesInSequence().get(i);
        QuantificationResult result = LDAResultReader.readResultFile(analysisModule_.getFullFilePath(exp).getAbsolutePath(),  new Hashtable<String,Boolean>());
        Hashtable<String,Hashtable<String,Vector<LipidomicsMSnSet>>> resultsOfExp = new Hashtable<String,Hashtable<String,Vector<LipidomicsMSnSet>>>();
        for (String className : result.getIdentifications().keySet()){
          Hashtable<String,Vector<LipidomicsMSnSet>> msnFound = new Hashtable<String,Vector<LipidomicsMSnSet>>();
          Vector<LipidParameterSet> sets = result.getIdentifications().get(className);
          for (LipidParameterSet set : sets){
            if (!(set instanceof LipidomicsMSnSet)) continue;
            Vector<LipidomicsMSnSet> ofOneAnalyte = new Vector<LipidomicsMSnSet>();
            if (msnFound.containsKey(set.getNameStringWithoutRt())) ofOneAnalyte = msnFound.get(set.getNameStringWithoutRt());
            ofOneAnalyte.add((LipidomicsMSnSet)set);
            msnFound.put(set.getNameStringWithoutRt(), ofOneAnalyte);
          }
          resultsOfExp.put(className, msnFound);
        }
        results.put(exp, resultsOfExp);
        headerLine += "\t\""+analysisModule_.getFullFilePath(exp).getName()+"\"";
      }
      
      stream.write((headerLine+"\n").getBytes());
      Hashtable<String,Vector<String>> correctAnalyteSequence = this.analysisModule_.getAllMoleculeNames();
      
      for (String className:this.heatmaps_.keySet()) {
        LinkedHashMap<String,String> mods = new LinkedHashMap<String,String>();
        for (String molName: correctAnalyteSequence.get(className)){
          //first key is the lipid molecular species identifier (FA level), second key is the experiment, third key is the modification
          LinkedHashMap<String,Hashtable<String,Hashtable<String,Vector<LipidomicsMSnSet>>>> sortedByMolSpecies = new  LinkedHashMap<String,Hashtable<String,Hashtable<String,Vector<LipidomicsMSnSet>>>>();
          for (int i=0; i!=analysisModule_.getExpNamesInSequence().size(); i++) {
            String exp = analysisModule_.getExpNamesInSequence().get(i);
            if (!results.containsKey(exp) || !results.get(exp).containsKey(className) || !results.get(exp).get(className).containsKey(molName)) continue;
            for (LipidomicsMSnSet set : results.get(exp).get(className).get(molName)){
              for (Object msnNames : set.getMSnIdentificationNames()){    
                String nameString = "";
                String faId = "";
                if (msnNames instanceof Vector){
                  for (String name : (Vector<String>)msnNames){
                    nameString += name+"|";
                  }
                  faId = nameString.substring(0,nameString.indexOf("|"));
                  nameString = nameString.substring(0,nameString.length()-1);
                }else{
                  nameString = (String)msnNames;
                  faId = nameString;
                }
                //this command has to be commented for position specific detection
                faId = faId.replaceAll("/", "_");
//                if (faId.indexOf("-_")>-1) faId = faId.replaceAll("-_", "");
//                else if (faId.indexOf("_-")>-1) faId = faId.replaceAll("_-", "");
                if (className.equalsIgnoreCase("DG") && set.getModificationName().equalsIgnoreCase("Na") && faId.split("_").length==2)
                  faId+="_-"; 
                Hashtable<String,Hashtable<String,Vector<LipidomicsMSnSet>>> oneMolSpecies = new Hashtable<String,Hashtable<String,Vector<LipidomicsMSnSet>>>();
                if (sortedByMolSpecies.containsKey(faId)) oneMolSpecies = sortedByMolSpecies.get(faId);
                else{
                  for (String otherId : sortedByMolSpecies.keySet()){
                    //TODO: this might not be correct - check when doing MAF-exports
                    if (StaticUtils.isAPermutedVersion(faId, otherId,LipidomicsConstants.CHAIN_SEPARATOR_NO_POS)){
                      faId = otherId;
                      oneMolSpecies = sortedByMolSpecies.get(faId);
                      break;
                    }
                  }
                }
                Hashtable<String,Vector<LipidomicsMSnSet>> speciesPerExp = new Hashtable<String,Vector<LipidomicsMSnSet>>();
                if (oneMolSpecies.containsKey(exp)) speciesPerExp = oneMolSpecies.get(exp);
                Vector<LipidomicsMSnSet> foundPerMod = new Vector<LipidomicsMSnSet>();
                if (speciesPerExp.containsKey(set.getModificationName())) foundPerMod = speciesPerExp.get(set.getModificationName());
                foundPerMod.add(set);
                speciesPerExp.put(set.getModificationName(), foundPerMod);
                if (!mods.containsKey(set.getModificationName())) mods.put(set.getModificationName(), set.getModificationName());
                oneMolSpecies.put(exp, speciesPerExp);
                sortedByMolSpecies.put(faId, oneMolSpecies);
              }
            } 
          }
          // if there are species present where no FAs are known - add them to the other ones
          if (sortedByMolSpecies.containsKey(molName) && sortedByMolSpecies.size()>1){
            Hashtable<String,Hashtable<String,Vector<LipidomicsMSnSet>>> totalMolSpecies = sortedByMolSpecies.get(molName);
            Hashtable<String,String> foundModsWithStructure = new Hashtable<String,String>();
            for (Hashtable<String,Hashtable<String,Vector<LipidomicsMSnSet>>> faResults : sortedByMolSpecies.values()){
              for (Hashtable<String,Vector<LipidomicsMSnSet>> expResults : faResults.values()){
                for (String mod : expResults.keySet()) foundModsWithStructure.put(mod, mod);
              }
            }
            Hashtable<String,String> foundMods = new Hashtable<String,String>();
            for (Hashtable<String,Vector<LipidomicsMSnSet>> expResults : totalMolSpecies.values()){
              for (String mod : expResults.keySet()) foundMods.put(mod, mod);
            }
            for (String faId : sortedByMolSpecies.keySet()){
              if (faId.equalsIgnoreCase(molName)) continue;
              Hashtable<String,Hashtable<String,Vector<LipidomicsMSnSet>>> oneMolSpecies = sortedByMolSpecies.get(faId);
              Hashtable<String,String> modsToAdd = new Hashtable<String,String>();
              for (Hashtable<String,Vector<LipidomicsMSnSet>> expResults : oneMolSpecies.values()){
                for (String mod : expResults.keySet()) modsToAdd.put(mod, mod);
              }
              for (String exp : totalMolSpecies.keySet()){
                Hashtable<String,Vector<LipidomicsMSnSet>> totalPerExp = totalMolSpecies.get(exp);
                if (oneMolSpecies.containsKey(exp)){
                  Hashtable<String,Vector<LipidomicsMSnSet>> speciesPerExp = oneMolSpecies.get(exp);
                  for (String mod : totalPerExp.keySet()){
                    if (foundModsWithStructure.containsKey(mod) && !modsToAdd.containsKey(mod)) continue;
                    Vector<LipidomicsMSnSet> totalPerMod = totalPerExp.get(mod);
                    if (speciesPerExp.containsKey(mod)){
                      Vector<LipidomicsMSnSet> foundPerMod = speciesPerExp.get(mod);
                      for (LipidomicsMSnSet hit : totalPerMod){
                        int count = 0;
                        for (LipidomicsMSnSet otherHit : foundPerMod){
                          if (Float.parseFloat(otherHit.getRt())<Float.parseFloat(hit.getRt()))count++;
                          else break;
                        }
                        foundPerMod.add(count,hit);
                      }
                      speciesPerExp.put(mod, foundPerMod);
                    }else{
                      speciesPerExp.put(mod, totalPerMod);
                    }
                  }
                  oneMolSpecies.put(exp,speciesPerExp);
                }else{
                  Hashtable<String,Vector<LipidomicsMSnSet>> toAdd = new Hashtable<String,Vector<LipidomicsMSnSet>>();
                  for (String mod : foundMods.keySet()){
                    if (!totalPerExp.containsKey(mod)) continue;
                    if (!foundModsWithStructure.containsKey(mod) || modsToAdd.containsKey(mod)) toAdd.put(mod, totalPerExp.get(mod));
                  }
                  if (toAdd.size()>0) oneMolSpecies.put(exp, toAdd);
                }
              }
              sortedByMolSpecies.put(faId,oneMolSpecies);
            }
            sortedByMolSpecies.remove(molName);
          }
          
          for (String faId : sortedByMolSpecies.keySet()){
            //System.out.println(className+faId);
            for (String mod : mods.keySet()){
              boolean modFound = false;
              String dbId = className+" "+molName;
              String chemicalFormula = "";
              String formattedMod = "";
              String metaboIdent = "";
              String fragmentString = "";
              String chargeString = "";
              String rtString = "";
              Hashtable<String,Double> totalMasses = new Hashtable<String,Double>();
              Hashtable<String,Integer> totalMassCount = new Hashtable<String,Integer>();
              Hashtable<String,Integer> posDefined = new Hashtable<String,Integer>();
              Hashtable<String,Double> abundances = new Hashtable<String,Double>();
              for (int i=0; i!=analysisModule_.getExpNamesInSequence().size(); i++){
                String exp  = analysisModule_.getExpNamesInSequence().get(i);
                Hashtable<String,Integer> posExpDefined = new Hashtable<String,Integer>();
                Hashtable<String,Integer> faDefined = new Hashtable<String,Integer>();
                Hashtable<String,Integer> nothingDefined = new Hashtable<String,Integer>();
                if (i!=0){
                  metaboIdent += "|";
                  fragmentString += "|";
                  rtString += "|";
                }
                double abundanceOfAnalyte = 0d;
                List<Float> fragments = new ArrayList<Float>();
                if (sortedByMolSpecies.get(faId).containsKey(exp) && sortedByMolSpecies.get(faId).get(exp).containsKey(mod)){
                  for (LipidomicsMSnSet set : sortedByMolSpecies.get(faId).get(exp).get(mod)){
                    for (Object msnNames : set.getMSnIdentificationNames()){    
                      String nameString = "";
                      String detectedFaId = "";
                      if (msnNames instanceof Vector){
                        for (String name : (Vector<String>)msnNames){
                          nameString += name+"|";
                        }
                        //the replaceAll has to be removed for position specific comparison
                        detectedFaId = nameString.substring(0,nameString.indexOf("|")).replaceAll("/", "_");
                        nameString = nameString.substring(0,nameString.length()-1);
                      }else{
                        nameString = (String)msnNames;
                      //the replaceAll has to be removed for position specific comparison
                        detectedFaId = nameString.replaceAll("/", "_");
                      }
                      boolean isCorrect = false;
                      if (className.equalsIgnoreCase("DG") && set.getModificationName().equalsIgnoreCase("Na") && detectedFaId.split("_").length==2){
                        detectedFaId+="_-"; 
                      }
                      //TODO: this might not be correct - check when doing MAF-exports
                      if (detectedFaId.equalsIgnoreCase(faId) || StaticUtils.isAPermutedVersion(detectedFaId, faId,LipidomicsConstants.CHAIN_SEPARATOR_NO_POS) || detectedFaId.equalsIgnoreCase(molName)) isCorrect = true;
                      if (!isCorrect)continue;
                      if (chemicalFormula.length()==0){
                        chemicalFormula = set.getAnalyteFormula().replaceAll(" ", "");
                        formattedMod = set.getModificationFormula().replaceAll(" ", "");
                        if (formattedMod.startsWith("-")) {
                          formattedMod = formattedMod.replaceAll("-", "");
                          formattedMod = "-" + formattedMod;
                        } else {
                          formattedMod = "+" + formattedMod;
                        }
                      }
                      //System.out.println("              "+nameString+" ; "+exp+" ; "+mod+" ; "+set.getRt());
                      modFound = true;
                      if (nameString.contains("/")){
                        int count = 0;
                        int count2 = 0;
                        if (posDefined.containsKey(nameString)) count = posDefined.get(nameString);
                        if (posExpDefined.containsKey(nameString)) count2 = posExpDefined.get(nameString);
                        count++;
                        count2++;
                        posDefined.put(nameString,count);
                        posExpDefined.put(nameString, count2);
                      } else if (nameString.contains("_")){
                        int count = 0;
                        if (nameString.contains("_") && faDefined.containsKey(nameString)) count=faDefined.get(nameString);
                        else if (nothingDefined.containsKey(nameString)) count = nothingDefined.get(nameString);
                        count++;
                        if (nameString.contains("_")) faDefined.put(nameString,count);
                        else nothingDefined.put(nameString,count);
                      }
                      double totalMass = 0d;
                      int totalCount = 0;
                      if (totalMasses.containsKey(mod)){
                        totalMass = totalMasses.get(mod);
                        totalCount = totalMassCount.get(mod);
                      }
                      for (CgProbe probe : set.getIsotopicProbes().get(0)){
                        totalMass += (double)probe.Mz;
                        totalCount++;
                      }
                      totalMasses.put(mod, totalMass);
                      totalMassCount.put(mod, totalCount);
                      
                      for (CgProbe probe : set.getHeadGroupFragments().values()){
                        fragments.add(probe.Mz);
                      }
                      //abundance is currently not used for the fragments
                      double abundance = 0d;
                      if (set.getStatus()>LipidomicsMSnSet.HEAD_GROUP_DETECTED){
                        String faName = nameString.replaceAll("/","_");
                        if (faName.indexOf("|")>-1) faName = faName.substring(0,faName.indexOf("|"));
                        //TODO: this might not be correct - check when doing MAF-exports
                        Vector<String> fas = StaticUtils.splitChainCombiToEncodedStrings(faName,LipidomicsConstants.CHAIN_SEPARATOR_NO_POS);
                        for (int j=0; j!= fas.size(); j++){
                          String fa = fas.get(j);
                          if (!set.getChainFragments().containsKey(fa)) continue;
                          for (CgProbe probe : set.getChainFragments().get(fa).values()){
                            fragments.add(probe.Mz);
                          }
                        }
                        abundance = ((double)set.getArea())*set.getRelativeIntensity(nameString);
                      } else {
                        abundance = (double)set.getArea();
                      }
                      if (rtString.length()>0 && !rtString.endsWith("|")) rtString+=",";
                      rtString += set.getRt();
                      abundanceOfAnalyte += abundance;
                    }
                    if (chargeString.length()==0) chargeString = set.getCharge().toString();
                    
                  }
                  if (posExpDefined.size()>0){
                    boolean isATie = true;
                    int highestFound = 0;
                    String bestId = "";
                    for (String pos : posExpDefined.keySet()){
                      int count = posExpDefined.get(pos);
                      if (count>highestFound){
                        bestId = pos;
                        isATie = false;
                        highestFound = count;
                      } else if (count==highestFound) isATie = true;
                    }
                    if (posExpDefined.size()==0 || isATie) metaboIdent += className+" "+faId;
                    else metaboIdent += className+" "+bestId.replaceAll("\\|", ",");
                  } else if (faDefined.size()>0) metaboIdent += className+" "+faId;
                  else metaboIdent += className+" "+molName;
                  
                  Collections.sort(fragments);
                  Hashtable<String,String> used = new Hashtable<String,String>();
                  for (Float fragment : fragments){
                    String mz = Calculator.FormatNumberToString(fragment,2);
                    if (used.contains(mz)) continue;
                    if (fragmentString.length()>0 && !fragmentString.endsWith("|")) fragmentString+=",";
                    fragmentString+=mz;
                    used.put(mz, mz);
                  }
                }else{
                  metaboIdent += "-";
                  fragmentString += "-";
                  rtString += "-";
                }
                if (abundanceOfAnalyte>0) abundances.put(exp, abundanceOfAnalyte);              
              }
              if (modFound){
                String line = "";
                if (!faId.equalsIgnoreCase(molName)){
                  dbId+=","+className+" ";
                  boolean isATie = true;
                  int highestFound = 0;
                  String bestId = "";
                  for (String pos : posDefined.keySet()){
                    int count = posDefined.get(pos);
                    if (count>highestFound){
                      bestId = pos;
                      isATie = false;
                      highestFound = count;
                    } else if (count==highestFound) isATie = true;
                  }
                  if (posDefined.size()==0 || isATie || bestId.indexOf("|")!=-1){
                    String toAdd = faId;
                    if (toAdd.indexOf("-_")!=-1) toAdd = toAdd.replaceAll("-_", "");
                    if (toAdd.indexOf("_-")!=-1) toAdd = toAdd.replaceAll("_-", "");
                    dbId += toAdd;
                  } else dbId += bestId;
                }
                line += "\""+dbId+"\"";
                line += "\t\""+chemicalFormula+"\"";
                line += "\t\"\"";
                line += "\t\"\"";
                line += "\t\""+metaboIdent+"\"";
                double mass = totalMasses.get(mod)/((double)totalMassCount.get(mod));
                line += "\t\""+Calculator.FormatNumberToString(mass,4)+"\"";
                line += "\t\""+fragmentString+"\"";
                line += "\t\""+"CHEMMOD:"+formattedMod+"\"";
                line += "\t\""+chargeString+"\"";
                line += "\t\""+rtString+"\"";
                if (LipidomicsConstants.getMzTabSample().getSpecies().get(0).getCvAccession()!=null && LipidomicsConstants.getMzTabSample().getSpecies().get(0).getCvAccession().length()>0 &&
                    LipidomicsConstants.getMzTabSample().getSpecies().get(0).getName()!=null && LipidomicsConstants.getMzTabSample().getSpecies().get(0).getName().length()>0){
                  line += "\t\""+LipidomicsConstants.getMzTabSample().getSpecies().get(0).getCvAccession()+"\"";
                  line += "\t\""+LipidomicsConstants.getMzTabSample().getSpecies().get(0).getName()+"\"";
                }else{
                  line += "\t\""+"\"";//+"taxid";
                  line += "\t\""+"\"";//+"species";                  
                }
                line += "\t\""+"\"";//+"database";
                line += "\t\""+"\"";//+"database_version";
                line += "\t\""+"\"";//+"reliability";
                line += "\t\""+"\"";//+"uri";
                line += "\t\""+"LipidDataAnalyzer, "+Settings.VERSION+"\"";
                line += "\t\""+"\"";//+"search_engine_score";
                line += "\t\""+Calculator.mean(new Vector<Double>(abundances.values()))+"\"";
                String stdev = "NaN";
                String stderr = "NaN";
                if (abundances.size()>1){
                  stdev = String.valueOf(Calculator.stddeviation(new Vector<Double>(abundances.values())));
                  stderr = String.valueOf(Calculator.stddeviation(new Vector<Double>(abundances.values()))/Math.sqrt(abundances.size()));
                }
                line += "\t\""+stdev+"\"";
                line += "\t\""+stderr+"\"";
                for (int i=0; i!=analysisModule_.getExpNamesInSequence().size(); i++) {
                  String exp = analysisModule_.getExpNamesInSequence().get(i);
                  line += "\t\"";
                  if (abundances.containsKey(exp)) line+= abundances.get(exp);
                  line += "\"";
                }
                stream.write((line+"\n").getBytes());
              }
            }
          }
            
        }
      }
      stream.close();
    }
    catch (IOException e) {
      e.printStackTrace();
    }
    catch (ExcelInputFileException e) {
      e.printStackTrace();
    }
    catch (LipidCombinameEncodingException e) {
      e.printStackTrace();
    }     
  }
  

  //TODO: for SILDA analysis, remove when done
  public void exportSummary(File exportFile, boolean isGrouped)
  {
  	Hashtable<String, HeatMapDrawing> heatMaps = isGrouped ? groupHeatmaps_ : heatmaps_;
  	OmegaCollector omegaCollector = new OmegaCollector(); //to collect class overarching data
  	try (BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(exportFile));)
  	{
  		XSSFWorkbook workbook = new XSSFWorkbook();
    	for (String molGroup : heatMaps.keySet())
    	{
    		Sheet sheet = workbook.createSheet(molGroup);
    		HeatMapDrawing drawing = heatMaps.get(molGroup);
    		drawing.exportSummary(omegaCollector, sheet, workbook, out);
    	}
    	ExcelAndTextExporter.writeHeatMapData(omegaCollector, workbook, out);
    	workbook.write(out);
  	}
  	catch (NumberFormatException e) {new WarningMessage(new JFrame(), "Error", e.getMessage());}
    catch (FileNotFoundException e) {new WarningMessage(new JFrame(), "Error", e.getMessage());}
    catch (IOException e) {new WarningMessage(new JFrame(), "Error", e.getMessage());}
  }
  
  
  private int getMaxProcessors(){
    int maxProcessors = Runtime.getRuntime().availableProcessors();
    if (maxProcessors<1) maxProcessors = Integer.MAX_VALUE;
    return maxProcessors;
  }
  
  private int getAmountOfProcessorsPreferred(){
    int maxProcessors = getMaxProcessors();
    if (maxProcessors == Integer.MAX_VALUE) maxProcessors=1;
    int procs = maxProcessors-1;
    if (procs<1) procs=1;
    return procs;
  }
  
  private Hashtable<String,LipidParameterSet> getParamByAnalyteName(String analyteName, Vector<LipidParameterSet> params, Hashtable<String,String> modHash){
    Hashtable<String,LipidParameterSet> paramOfInterest = new Hashtable<String,LipidParameterSet>();
    String[] nameAndRt = StaticUtils.extractMoleculeRtAndModFromMoleculeName(analyteName);
    if (nameAndRt[1]==null) nameAndRt[1] = "";
    String[] nameAndRtParam1 = StaticUtils.extractMoleculeRtAndModFromMoleculeName(analyteName);
    if (nameAndRtParam1[1]==null) nameAndRtParam1[1]="";
    for (LipidParameterSet param : params){
      String[] nameAndRtParam2 = StaticUtils.extractMoleculeRtAndModFromMoleculeName(param.getNameString());
      if (nameAndRtParam2[1]==null) nameAndRtParam2[1]="";
      if (nameAndRtParam1[0].equalsIgnoreCase(nameAndRtParam2[0])){
        if (analysisModule_.neglectRtInformation(analyteName) || analysisModule_.isWithinRtGroupingBoundaries(Double.valueOf(nameAndRtParam1[1]), Double.valueOf(nameAndRtParam2[1]))){
          if (modHash.containsKey(param.getModificationName())){
            paramOfInterest.put(param.getModificationName(), param);
          }
        }  
      }  
    }
    return paramOfInterest;
  }
  
  /**
   * Cleans up fields related to results loaded by the statistical analysis.
   */
  private void cleanupResultView(){
  	LipidParameterSet.setOmegaInformationAvailable(false); //this ensures the export settings only give the option to export omega positions when any are present
  	groupDisplayNamesLookup_ = null;
  	expDisplayNamesLookup_ = null;
  	if (heatmaps_!=null){
      for (HeatMapDrawing map : heatmaps_.values()) map.cleanup();
      heatmaps_ = null;
    }
  	if (groupHeatmaps_!=null){
      for (HeatMapDrawing map : groupHeatmaps_.values()) map.cleanup();
      groupHeatmaps_ = null;
    }
  	molBarCharts_ = null;
  	if (colorChooserDialog_ != null) colorChooserDialog_.cleanup();
  	colorChooserDialog_ = null;
  	exportSettings_ = null;
  	exportSettingsGroup_ = null;
  	if (classOverviewPanel_ != null) classOverviewPanel_.cleanup();
  	classOverviewPanel_ = null;
  	if (classOverviewGroupPanel_ != null) classOverviewGroupPanel_.cleanup();
  	classOverviewGroupPanel_ = null;
    removeResultTabComponentsExceptFirst();
  }
  
  /**
   * changes the list sort type of the "Display results" table
   */
  public void changeListSorting(int sortType){
    orderResultsType_.put(this.currentSelectedSheet_, sortType);
    updateResultListSelectionTable();
    this.currentSelected_ = -1;
  }
  
  
  /**
   * returns start and stop times for MSn display in the 3D viewer
   * @param retTimes the retention time look up
   * @param start start scan of the found hits in the chromatogram
   * @param stop stop scan of the found hits in the chromatogram
   * @param probes the found MS1 probes to be displayed
   * @return float[0] = start retention time; float[1] = stop retention time
   */
  private float[] getMSn3DTimeStartStopValues(Hashtable<Integer,Float> retTimes, int start, int stop, Vector<CgProbe> probes){
    float[] startStop = new float[2];
    startStop[0] = retTimes.get(start);
    startStop[1] = retTimes.get(stop);
    if (probes!=null && probes.size()>0){
      float lowestTime = Float.MAX_VALUE;
      float highestTime = 0f;
      for (CgProbe probe : probes){
        if (probe.LowerValley<lowestTime) lowestTime = probe.LowerValley;
        if (probe.UpperValley>highestTime) highestTime = probe.UpperValley;
      }
      lowestTime -= Settings.getMsn3DDisplayTolerance();
      highestTime += Settings.getMsn3DDisplayTolerance();
      if (lowestTime<startStop[1] && lowestTime>startStop[0]) startStop[0] = lowestTime;
      if (highestTime>startStop[0]&& highestTime<startStop[1]) startStop[1] = highestTime;
    }
    return startStop;
  }
  
  private LipidParameterSet refreshSpectrumPainterRDI() throws RulesException, NoRuleException, IOException, SpectrummillParserException, CgException, HydroxylationEncodingException, ChemicalFormulaException, LipidCombinameEncodingException{
    LipidParameterSet param = msnUserInterfaceObject_.testForMSnDetection(getMs2LevelSpectrumSelected());
    Hashtable<Integer,Vector<RangeColor>> rangeColors = StaticUtils.createRangeColorVOs(param, ((LipidomicsTableModel)displayTable_.getModel()).getMSnIdentificationName(currentSelected_),
        result_.getFaHydroxyEncoding(), result_.getLcbHydroxyEncoding(), areTheseAlex123MsnFragments());
    if (rangeColors!=null) spectrumPainter_.refresh(param,rangeColors);
    else spectrumPainter_.clearRangeColors();
    this.spectrumSelected_.setText(spectrumPainter_.getSpectSelectedText());
    this.rtSelected_.setText(spectrumPainter_.getRtSelectedText());
    displayPrecursorMasses(spectrumPainter_.getPrecursorMassSelected());
    this.msLevelSelected_.setText(spectrumPainter_.getMsLevelSelected());
    return param;
  }
  
  /**
   * Paints a new spectra with certain parameters
   * @param params
   * @throws CgException
   */
  public void updateSpectra(LipidParameterSet params, boolean newPrec) throws CgException
  {
    displayTolerancePanel_.getLockMzRange().setSelected(false);
    if(newPrec ==  true){
      //PAINTS THE NEW SPEKTRA WITH THE LPS
      float m2dGain = spectrumPainter_.getM2dGain();   
      spectrumPanel_.remove(spectrumPainter_);   
      //int msLevel = 2;
      float tol = LipidomicsConstants.getMs2PrecursorTolerance(params.Mz[0]);
      Hashtable<Integer,Vector<String>> spectraRaw = reader_.getMsMsSpectra(params.Mz[0]-tol, params.Mz[0]+tol, -1f, -1f);
      int[] borders = null;
      for (Integer level : spectraRaw.keySet()){
        //TODO: a dummy call to set the values for getBordersOfLastMs2Extraction()
        reader_.translateSpectraToChroms(spectraRaw.get(level), level, LipidomicsConstants.getMs2ChromMultiplicationFactorForInt());
        int[] bordersLevel = reader_.getBordersOfLastMs2Extraction();
        if (borders==null){
          borders = bordersLevel;
        }else{
          if (bordersLevel[0]<borders[0]) borders[0] = bordersLevel[0];
          if (bordersLevel[2]<borders[2]) borders[2] = bordersLevel[2];
          if (bordersLevel[1]>borders[1]) borders[1] = bordersLevel[1];
          if (bordersLevel[3]>borders[3]) borders[3] = bordersLevel[3];
        }
      }
      if (borders==null){
        throw new CgException("There are no spectra present!");
      }      

      
      @SuppressWarnings("rawtypes")
      Vector<Hashtable> rtNrSpectraAndPrecursor = reader_.getRtNrSpectrumHash(spectraRaw);
      @SuppressWarnings("unchecked")
      Hashtable<Integer,String> rtNrSpectrumHash = (Hashtable<Integer,String>)rtNrSpectraAndPrecursor.get(0);
      @SuppressWarnings("unchecked")
      Hashtable<Integer,Vector<Double>> rtNrPrecursorHash = (Hashtable<Integer,Vector<Double>>)rtNrSpectraAndPrecursor.get(1);
      @SuppressWarnings("unchecked")
      Hashtable<Integer,Integer> scanNrLevelHash = (Hashtable<Integer,Integer>)rtNrSpectraAndPrecursor.get(2);
      Hashtable<Integer,Float> retTimes = reader_.getMsmsRetentionTimes();
      float peakRt = Float.parseFloat(params.getRt())*60f;

      int extendedStart = (borders[0]*9)/10;
      int extendedStop = (borders[1]*21)/20;
      Lipidomics2DSpectraChromPainter spectrumPainter = new Lipidomics2DSpectraChromPainter(analyzer_,rtNrSpectrumHash, rtNrPrecursorHash, scanNrLevelHash, retTimes,
          peakRt, extendedStart,extendedStop,LipidomicsConstants.getMs2ChromMultiplicationFactorForInt(),tol*2,this,
          extendedStart,extendedStop, true,new Vector<CgProbe>(),new Vector<CgProbe>(),1,1, this.relAbund_.isSelected(),
          params,new Hashtable<Integer,Vector<RangeColor>>(),new Double(annotationThreshold_.getText()),false);      

      //SETS THE ZOOM BACK BEFORE THE REFRESHMENT
      if(mz_minTimeText_.getText() != null){
        if(!(mz_minTimeText_.getText().isEmpty())){
          spectrumPainter.setMinDispTime2d(Float.parseFloat(mz_minTimeText_.getText()));
        }
      }
      if(mz_maxTimeText_.getText() != null){
        if(!(mz_maxTimeText_.getText().isEmpty())){
          spectrumPainter.setMaxDispTime2d(Float.parseFloat(mz_maxTimeText_.getText()));
        }
      }
  
      //SETS THE GAIN BACK BEFORE THE REFRESHMENT
      spectrumPainter.setM2dGain(m2dGain);     
      spectrumPainter_ = spectrumPainter;     
    
      spectrumPanel_.add(spectrumPainter_,BorderLayout.CENTER); 
      spectrumPainter.setBackground(Color.WHITE);       
    
      displayPrecursorMasses(spectrumPainter_.getPrecursorMassSelected());
      rtSelected_.setText(spectrumPainter_.getRtSelectedText());
      msLevelSelected_.setText(spectrumPainter_.getMsLevelSelected());
      spectrumPainter_.draw2DDiagram(extendedStart,extendedStop, true);  
    
      LipidParameterSet pr = params;
      try {
        pr = refreshSpectrumPainterRDI();
      } catch (RulesException | NoRuleException | IOException
        | SpectrummillParserException | HydroxylationEncodingException | ChemicalFormulaException | LipidCombinameEncodingException e1) {
      }
    
      //REFRESHES THE CHROM PAINTER FROM HERE ONLY FOR NEW PRECUSORS
      int charge = 1;
      if (pr.ProbeCount()>0)
        charge = pr.Probe(0).Charge;
      else if (pr.getCharge()!=null&&params_.getCharge()>1)
        charge = pr.getCharge();
      float currentIsotopicMass = pr.Mz[0];    
      isotope_.setSelectedIndex(0);
      float startFloat = currentIsotopicMass-Float.parseFloat(this.displayTolerancePanel_.getDisplayMinusTolerance().getText());
      float stopFloat = currentIsotopicMass+Float.parseFloat(this.displayTolerancePanel_.getDisplayPlusTolerance().getText());
      Vector<CgProbe> storedProbes = new Vector<CgProbe>();
      Vector<CgProbe> selectedProbes = new Vector<CgProbe>();
      //Vector<CgProbe> storedProbes = ms1ProbesWhileMs2Display_.get(0);
      //Vector<CgProbe> selectedProbes = ms1ProbesWhileMs2Display_.get(1);
      l2DPainter_.getGraphics().dispose();
      l2dPanel_.remove(l2DPainter_);
      l2DPainter_ = null;
  
      try {
        String[] rawLines = reader_.getRawLines(startFloat, stopFloat, result_.getMsLevels().get(currentSelectedSheet_));
        Hashtable<Integer,Float> rtTimes = reader_.getRetentionTimes();
        Lipidomics2DPainter l2DPainter = new Lipidomics2DPainter(analyzer_,rawLines, rtTimes,
            startFloat,stopFloat,reader_.getMultiplicationFactorForInt_()/reader_.getLowestResolution_(),pr.LowerMzBand*2,this,
            MSMapViewer.DISPLAY_TIME_MINUTES,currentIsotopicMass-pr.LowerMzBand, currentIsotopicMass+pr.UpperMzBand, false,
            storedProbes,selectedProbes,Integer.parseInt((String)this.isotope_.getSelectedItem()),charge,result_.getMsLevels().get(currentSelectedSheet_),
            shotgunIsDisplayed_);
        l2DPainter.preChromatogramExtraxtion(currentIsotopicMass-pr.LowerMzBand, currentIsotopicMass+pr.UpperMzBand);
        l2DPainter_ = l2DPainter;
        l2DPainter_.setBackground(Color.WHITE);
        l2dPanel_.add(l2DPainter,BorderLayout.CENTER);
      } catch (CgException e) {
        e.printStackTrace();
      }
      l2DPainter_.repaint();
      l2DPainter_.paintMs2Position(Float.parseFloat(rtSelected_.getText()));
      majorSplitPane_.repaint(); 
      topSplitPane_.setResizeWeight(0.58);
    
      showMs2InputElements();
      l2dPanel_.invalidate();
      l2dPanel_.updateUI();
      spectrumPanel_.invalidate();
      spectrumPanel_.updateUI();
      
    } else {
      try {
        refreshSpectrumPainterRDI();
      }
      catch (RulesException | NoRuleException | IOException
          | SpectrummillParserException | HydroxylationEncodingException | ChemicalFormulaException | LipidCombinameEncodingException e) {
        e.printStackTrace();
      }
    }
  }
  
  /**
   * Returns the selected spectra
   * @return
   */
  public int getSelectedSpectrumNumber()
  {
    return spectrumPainter_.getSpectrumSelected();
  }
  
  /**
   * returns the current number of the MS2 spectrum from a the list of sequentially ordered MS2-spectra 
   */
  public int getMs2LevelSpectrumSelected()
  {
    return spectrumPainter_.getMs2LevelSpectrumSelected();
  }

  
  /**
   * Creates the interface via the RuleDefinitionInterface class
   */
  public void newRule(int position) 
  {     
    LipidParameterSet data = getAnalyteInTableAtPosition(position);
    if (data!=null) {
      if (data instanceof LipidomicsMSnSet)
        try {
          data = new LipidomicsMSnSet((LipidomicsMSnSet)data);
        } catch (LipidCombinameEncodingException e) {
          e.printStackTrace();
        }
      else
        data = new LipidParameterSet(data);
    }
    currentSelected_ = position;
    try 
    {      
      MSnAnalyzer analyzer = new MSnAnalyzer(null,(String)selectedSheet_.getSelectedItem(),data.getModificationName(),data,analyzer_);      
      data = analyzer.getResult();
    }
    catch (RulesException e) {
      e.printStackTrace();
    }
    catch (IOException e) {
      e.printStackTrace();
    }
    catch (SpectrummillParserException e) {
      e.printStackTrace();
    }
    catch (CgException e) {
      e.printStackTrace();
    }
    catch (HydroxylationEncodingException | ChemicalFormulaException | LipidCombinameEncodingException e) {
      e.printStackTrace();
    }
    
    if (!showMs2(data,true)) return; 
    try {
      userInterface = new RuleDefinitionInterface((String)selectedSheet_.getSelectedItem(), data, analyzer_, reader_.getHighestMsLevel(), this);  
      if(viewer_!=null)
      {
        topSplitPane_.remove(viewer_);  
      } 
      topSplitPane_.setTopComponent(userInterface);
      majorSplitPane_.setDividerLocation(0.72);
      majorSplitPane_.setBottomComponent(this.spectrumPanel_);
      finalButtonSection_ = userInterface.makeFinalButtonSection();
      spectrumPanel_.add(finalButtonSection_,BorderLayout.SOUTH);
      msnUserInterfaceObject_ = userInterface;
      refreshSpectrumPainterRDI();
    }
    catch (NoRuleException | HydroxylationEncodingException | LipidCombinameEncodingException e){
      
    }
    catch (RulesException | IOException
        | SpectrummillParserException | CgException
        | ChemicalFormulaException e) {
      e.printStackTrace();
    }
  }
  
  public void updateCurrentSpectrum(){
    String cutoff = annotationThreshold_.getText();
    if (cutoff==null || cutoff.length()<1) return;
    try{
      double cutoffValue = Double.parseDouble(cutoff);
      spectrumPainter_.setAnnotationThreshold(cutoffValue);
    } catch (NumberFormatException nfx){}
  }
  
  private int[] getScrollPaneWidthAndHeight(HeatMapDrawing drawing){
    int[] wh = new int[2];
    wh[0] = drawing.getTotalSize().width+30;
    if (wh[0]<400) wh[0] = 360; 
    wh[1] = 700;
    if (this.getWidth()<(wh[0]+15)) wh[0] = this.getWidth()-15;
    if (this.getHeight()<(wh[1]+15)) wh[1] = this.getHeight()-15;
    return wh;
  }
  
  /**
   * @return does this class contain fragments originating from Alex123 MSn target lists
   */
  private boolean areTheseAlex123MsnFragments(){
    String className = (String)selectedSheet_.getSelectedItem();
    if (result_!=null && result_.getConstants()!=null && result_.getConstants().getAlexTargetlistUsed()!=null && result_.getConstants().getAlexTargetlistUsed().containsKey(className) && result_.getConstants().getAlexTargetlistUsed().get(className))
      return true;
    else
      return false;
  }
  
  /**
   * displays the precursor masses of a spectrum
   * @param precMasses the precursor masses
   */
  private void displayPrecursorMasses(Vector<Double> precMasses){
    for (JLabel label : precursorSelected_) label.setVisible(false);
    if (precMasses==null || precMasses.size()==0){
      precursorSelected_[0].setText("");
      precursorSelected_[0].setVisible(true);
      return;
    }
    int count = 0;
    for (int i=(precMasses.size()-1); i!=-1; i--){
      precursorSelected_[count].setText(String.valueOf(Calculator.roundDBL(precMasses.get(i),4)));
      precursorSelected_[count].setVisible(true);
      count++;
    }
  }
  
//  private void writeDisplayDataToExcelFormat(String[] rawLines,  Hashtable<Integer,Float> retentionTimes, float startFloat, float stopFloat){
//    String excelFile = "D:\\lipidomics\\20170925_Glucose\\rawData.xlsx";
//    float startTime = 6.2f*60f;
//    float stopTime = 7.9f*60f;
////    float startTime = 0f;
////    float stopTime = 0f;
//    
//    int startScan = -1;
//    int stopScan = 100000000;
//    if (startTime>0){
//      Enumeration<Integer> keys = retentionTimes.keys();
//      int lowestScan = 100000000;
//      while (keys.hasMoreElements()){
//        Integer key = keys.nextElement();
//        if (retentionTimes.get(key)>=startTime&&key<lowestScan)
//          lowestScan = key;
//      }
//      startScan = lowestScan;
//    }
//    if (stopTime>0){
//      Enumeration<Integer> keys = retentionTimes.keys();
//      int highestScan = 0;
//      while (keys.hasMoreElements()){
//        Integer key = keys.nextElement();
//        if (retentionTimes.get(key)<=stopTime&&key>highestScan)
//          highestScan = key;
//      }
//      stopScan = highestScan;
//    }
//
//    int dataAxisSize = retentionTimes.size();
//    if (stopScan<100000000){
//      if (dataAxisSize>stopScan)
//        dataAxisSize = stopScan+1;
//    }
//    if (startScan>-1){
//      dataAxisSize = dataAxisSize-startScan;
//    }
//
//    
//    float[][]data = new float[rawLines.length][dataAxisSize];
//    for (int i = 0; i != rawLines.length; i++) {
//      if (rawLines[i] != null && rawLines[i].length() > 0) {
//        ByteBuffer buffer = ByteBuffer.wrap(Base64.decode(rawLines[i]));
//        while (buffer.hasRemaining()) {
//          int scanNumber = buffer.getInt();
//          float intensity = buffer.getFloat();
//          if (scanNumber>=startScan&&scanNumber<=stopScan){
//            if (startScan>-1)
//              data[i][scanNumber-startScan] += intensity;
//            else
//              data[i][scanNumber] += intensity;
//          }
//        }
//      }
//    }
//    
//    try{
//      Workbook workbook = new XSSFWorkbook();
//      CellStyle headerStyle = workbook.createCellStyle();
//      org.apache.poi.ss.usermodel.Font arial12font = workbook.createFont();
//      arial12font.setBoldweight(org.apache.poi.ss.usermodel.Font.BOLDWEIGHT_BOLD);
//      arial12font.setFontName("Arial");
//      arial12font.setFontHeightInPoints((short)12);
//      headerStyle.setFont(arial12font);
//      headerStyle.setAlignment(CellStyle.ALIGN_CENTER);
//
//      BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(excelFile));
//      Sheet sheet = workbook.createSheet("Data");
//      int add = 0;
//      if (startScan>-1) add = startScan;
//      int count = 0;
//      Row row = sheet.createRow(count);
//      Cell cell;
//      count++;
//      
//      for (int j=0; j!=dataAxisSize; j++){
//        String rt = String.valueOf(Calculator.roundFloat(retentionTimes.get(j+add)/60f,3));
//        cell = row.createCell(j+1);
//        cell.setCellStyle(headerStyle);
//        cell.setCellValue(rt);
//      }
//      
//      for (int i=0; i!=data.length; i++){
//        String mzString = String.valueOf(Calculator.roundDBL(((double)startFloat)+i*0.0001d,4))+"-"+String.valueOf(Calculator.roundDBL(((double)startFloat)+(i+1)*0.0001d,4));
//        row = sheet.createRow(count);
//        count++;
//        cell = row.createCell(0);
//        cell.setCellStyle(headerStyle);
//        cell.setCellValue(mzString);
//        sheet.setColumnWidth(0, 6000);
//        
//        for (int j=0; j!=dataAxisSize; j++){
//          cell = row.createCell(j+1);
//          cell.setCellStyle(headerStyle);
//          cell.setCellValue(data[i][j]);
//        }
//       
//      }
//      
//      workbook.write(out);
//      out.close();
//      workbook.close();
//
//    }catch(Exception ex){
//      ex.printStackTrace();
//    }
//    
//  }
  
  /**
   * generates the quantification pair VOs for the progress display table of the quantification
   * this method respects whether there is polarity switched data and assigns the quantification files correspondingly
   * @param rawFiles the MS data files
   * @param quantFiles the quantification files
   * @return the MS data/quantification file pairs
   */
  public static Vector<RawQuantificationPairVO> generateQuantificationPairVOs(Vector<File> rawFiles, Vector<File> quantFiles){
    Vector<RawQuantificationPairVO> pairs = new Vector<RawQuantificationPairVO>();
    Vector<File> positivePolarity = new Vector<File>();
    Vector<File> negativePolarity = new Vector<File>();
    Hashtable<String,String> usedFiles = new Hashtable<String,String>();
    //extract the polarity switched files
    for (File rawFile: rawFiles){
      if (!rawFile.getAbsolutePath().endsWith(".chrom") || !rawFile.isDirectory())
        continue;
      File headerFile = new File(StringUtils.getChromFilePaths(rawFile.getAbsolutePath())[1]);
      if (headerFile.exists() && headerFile.isFile()){
        ChromatogramHeaderFileReader headerReader = new ChromatogramHeaderFileReader(headerFile.getAbsolutePath());
        try {
          headerReader.readFile();
          String polarity = headerReader.getPolaritySwitched();
          if (polarity==null) continue;
          if (polarity.equalsIgnoreCase(GlobalConstants.CHROMATOGRAM_HEADER_FILE_POLARITY_POSITIVE) ||
              polarity.equalsIgnoreCase(GlobalConstants.CHROMATOGRAM_HEADER_FILE_POLARITY_NEGATIVE)){
            usedFiles.put(rawFile.getAbsolutePath(), rawFile.getAbsolutePath());
            if (polarity.equalsIgnoreCase(GlobalConstants.CHROMATOGRAM_HEADER_FILE_POLARITY_POSITIVE))
              positivePolarity.add(rawFile);
            else if (polarity.equalsIgnoreCase(GlobalConstants.CHROMATOGRAM_HEADER_FILE_POLARITY_NEGATIVE))
              negativePolarity.add(rawFile);
          }
        } catch (IndexFileException | CgException e) {e.printStackTrace();}
      }
    }
    
    //exclude raw files that produced polarity switched files
    for (File rawFile: rawFiles){
      if (usedFiles.containsKey(rawFile.getAbsolutePath()))
        continue;
      String fileSuffix = rawFile.getAbsolutePath().substring(rawFile.getAbsolutePath().lastIndexOf("."));
      for (File posFile : positivePolarity){
        String comparePath = posFile.getAbsolutePath().substring(0,posFile.getAbsolutePath().length()-".chrom".length()-RawToChromThread.FILE_SUFFIX_POLARITY_POSITIVE.length());
        if (rawFile.getAbsolutePath().equalsIgnoreCase(comparePath+fileSuffix)){
          usedFiles.put(rawFile.getAbsolutePath(), rawFile.getAbsolutePath());
          break;
        }
      }
      for (File negFile : negativePolarity){
        String comparePath = negFile.getAbsolutePath().substring(0,negFile.getAbsolutePath().length()-".chrom".length()-RawToChromThread.FILE_SUFFIX_POLARITY_NEGATIVE.length());
        if (rawFile.getAbsolutePath().equalsIgnoreCase(comparePath+fileSuffix)){
          usedFiles.put(rawFile.getAbsolutePath(), rawFile.getAbsolutePath());
          break;
        }
      }

    }
    
    //add the positive ones
    for (File posFile : positivePolarity){
      boolean foundPositive = false;
      for (File quantFile: quantFiles){
        if (quantFile.getName().indexOf(GlobalConstants.CHROMATOGRAM_HEADER_FILE_POLARITY_POSITIVE)==-1)
          continue;
        pairs.add(new RawQuantificationPairVO(posFile,quantFile));
        foundPositive = true;
      }
      if (!foundPositive)
        new WarningMessage(new JFrame(), "Warning", "The file: "+posFile+" is of positive polarity, but no adequate quant file exist (must contain "+GlobalConstants.CHROMATOGRAM_HEADER_FILE_POLARITY_POSITIVE+" in the name)");
    }
    //add the negative ones
    for (File negFile : negativePolarity){
      boolean foundNegative = false;
      for (File quantFile: quantFiles){
        if (quantFile.getName().indexOf(GlobalConstants.CHROMATOGRAM_HEADER_FILE_POLARITY_NEGATIVE)==-1)
          continue;
        pairs.add(new RawQuantificationPairVO(negFile,quantFile));
        foundNegative = true;
      }
      if (!foundNegative)
        new WarningMessage(new JFrame(), "Warning", "The file: "+negFile+" is of negative polarity, but no adequate quant file exist (must contain "+GlobalConstants.CHROMATOGRAM_HEADER_FILE_POLARITY_NEGATIVE+" in the name)");
    }
    
    //add the rest
    for (File rawFile: rawFiles){
      if (usedFiles.containsKey(rawFile.getAbsolutePath()))
        continue;
      for (File quantFile: quantFiles){
        pairs.add(new RawQuantificationPairVO(rawFile,quantFile));
      }
    }
    return pairs;
  }

  public void recalculateMSn(int position)
  {
    LipidParameterSet data = new LipidParameterSet(getAnalyteInTableAtPosition(position));
    currentSelected_ = position;
    try {
      String className = (String)selectedSheet_.getSelectedItem();
      MSnAnalyzer analyzer = new MSnAnalyzer(null,className,data.getModificationName(),data,analyzer_);
      recalcDialog_ = new RecalculateMSnDialog(className,analyzer,this);
    } catch (RulesException | IOException | SpectrummillParserException | CgException | HydroxylationEncodingException | ChemicalFormulaException | LipidCombinameEncodingException e) {
      e.printStackTrace();
      new WarningMessage(new JFrame(), "ERROR", e.getMessage());
    }
  }

  public void editRt(int position){
    currentSelected_ = position;
    String rt = getAnalyteInTableAtPosition(position).getRt();
    editRtDialog_ = new EditRtDialog(rt,this);
  }
  
  public void editOmegaAssignment(int position) {
    currentSelected_ = position;
    String molecularSpecies = resultPositionToMolecularSpeciesLookup_.get(position);
    LipidParameterSet param = getAnalyteInTableAtPosition(position);
    try {
      if (Integer.parseInt(RulesContainer.getAmountOfChains(StaticUtils.getRuleName(currentSelectedSheet_, param.getModificationName()))) < 2) {
        molecularSpecies = param.getNameStringWithoutRt();
      } 
    } catch (RulesException | NoRuleException | IOException | SpectrummillParserException ex) {
      System.out.println(ex.getMessage());
    }
    
    Vector<FattyAcidVO> chainCombination = new Vector<FattyAcidVO>();
    if (molecularSpecies == null) {
      if (param instanceof LipidomicsMSnSet) {
        LipidomicsMSnSet lipidomicsMSnSet = (LipidomicsMSnSet)param;
        Set<String> unlabeledAssignmentsSet = lipidomicsMSnSet.getHumanReadableNameSet();
        if (!unlabeledAssignmentsSet.isEmpty()) {
          Object[] molecularSpeciesArray = unlabeledAssignmentsSet.toArray();
          if (unlabeledAssignmentsSet.size() == 1) {
            molecularSpecies = (String)molecularSpeciesArray[0];
          } else {
            JFrame frame = new JFrame();
            frame.setLocationRelativeTo(frame_);
            molecularSpecies = (String) JOptionPane.showInputDialog(frame, 
                "Select a molecular species!",
                "\u03C9-C=C assignment",
                JOptionPane.QUESTION_MESSAGE, 
                null, 
                molecularSpeciesArray, 
                molecularSpeciesArray[0]);
          }
        } 
      } else {
        new WarningMessage(new JFrame(), "Error", "This analyte lacks required MSn information for \u03C9 - double bond assignment.");
      } 
    } 
    if (molecularSpecies != null) {
      try {
      	//TODO: the chain mod separator should be known here => give it lipidomicsConstants? change the staticutils method?
        chainCombination = StaticUtils.decodeFAsFromHumanReadableName(molecularSpecies,result_.getFaHydroxyEncoding(), result_.getLcbHydroxyEncoding(), areTheseAlex123MsnFragments(),null);
      } catch (LipidCombinameEncodingException ex) { 
        ex.printStackTrace();
      }
      boolean doubleBondPresent = false;
      for (FattyAcidVO fattyAcid : chainCombination) {
        if (fattyAcid.getDoubleBonds() > 0) {
          doubleBondPresent = true;
        }
      }
      if (doubleBondPresent) {
        new EditOmegaAssignmentJTable(param, molecularSpecies, chainCombination, this, frame_);
      } else {
        new WarningMessage(new JFrame(), "Error", "This molecular species does not have any double bonds!");
      }
    } 
  }
  
  //TODO: rename if it is fine
  public void editOmegaAssignmentOld(int position) {
    currentSelected_ = position;
    LipidParameterSet param = getAnalyteInTableAtPosition(position);
    String rowToMolecularSpecies = resultPositionToMolecularSpeciesLookup_.get(position);
    Vector<FattyAcidVO> chainCombination = new Vector<FattyAcidVO>();
    if (rowToMolecularSpecies == null) {
      Set<String> molecularSpeciesSet = StaticUtils.getMolecularSpeciesSet(param.getOmegaInformation());
      if (molecularSpeciesSet.size() == 1) {
        rowToMolecularSpecies = molecularSpeciesSet.toString().replace("[", "").replace("]", "");
      } else if (!displayTolerancePanel_.getShowMSnNames().isSelected()){
        new WarningMessage(new JFrame(), "Warning", "Multiple MSn assignments for this analyte may exist. Check the box 'Show MSn' to edit the \u03C9 - double bond assignment for each.");
      }
    }
    if (rowToMolecularSpecies != null) {
      try {
      	//TODO: the chain mod separator should be known here => give it lipidomicsConstants? change the staticutils method?
        chainCombination = StaticUtils.decodeFAsFromHumanReadableName(rowToMolecularSpecies,result_.getFaHydroxyEncoding(), result_.getLcbHydroxyEncoding(), areTheseAlex123MsnFragments(),null);
      } catch (LipidCombinameEncodingException ex) { 
        ex.printStackTrace();
      }
      boolean doubleBondPresent = false;
      for (FattyAcidVO fattyAcid : chainCombination) {
        if (fattyAcid.getDoubleBonds() > 0) {
          doubleBondPresent = true;
        }
      }
      if (doubleBondPresent) {
        new EditOmegaAssignmentJTable(param, rowToMolecularSpecies, chainCombination, this, frame_);
      } else {
        new WarningMessage(new JFrame(), "Error", "This molecular species does not have any double bonds!");
      }
    }
  }
  
  /**
   * disables input fields specific for chromatography based data
   */
  private void disableChromatographyFeatures(){
    switchChromatographyFeatures(false);
  }

  /**
   * enables input fields specific for chromatography based data
   */
  private void enableChromatographyFeatures(){
    switchChromatographyFeatures(true);    
  }
  
  /**
   * copy content to clipboard
   */
  public void copyToClipboard() {
	  StringBuffer strBuf = new StringBuffer();
	  int sheets = selectedSheet_.getItemCount();
	  int selection = selectedSheet_.getSelectedIndex();
	
		ArrayList<String> al = new ArrayList<String>();
		for(int i = 0; i < sheets; i++) {
			selectedSheet_.setSelectedIndex(i);
			updateResultListSelectionTable();
			al.add((String) selectedSheet_.getSelectedItem());
		}
		Collections.sort(al);
		
	  for (int h = 0; h < al.size(); h++) {
		  
			selectedSheet_.setSelectedItem(al.get(h));
			updateResultListSelectionTable();
		  
		  int rowCount = displayTable_.getRowCount();
		  int columnCount = displayTable_.getColumnCount();
		  
		  strBuf.append(al.get(h) + "\n");
		  
		  for (int i = 0; i < rowCount; i++){
			for (int j = 0; j < columnCount; j++) {
				strBuf.append(displayTable_.getValueAt(i, j));
				strBuf.append("\t");
			}
		    
			strBuf.append("\n");
		  }
		  strBuf.append("\n");
	  }
	  selectedSheet_.setSelectedIndex(selection);
		
	  Toolkit toolkit = Toolkit.getDefaultToolkit();
	  Clipboard clipboard = toolkit.getSystemClipboard();
	  StringSelection strSel = new StringSelection(strBuf.toString());
	  clipboard.setContents(strSel, null);
	}
  
  /**
   * enables/disables input fields specific for chromatography based data
   * @param enable enabling when true, otherwise disabling
   */
  private void switchChromatographyFeatures(boolean enable){
    if (!enable){
      m_chkRaw_.setSelected(true);
    }
    m_chkSmooth_.setEnabled(enable);
    shotgunIsDisplayed_ = !enable;
  }
  
  public LinkedHashMap<String,Vector<String>> getSamplesOfGroups(){
    LinkedHashMap<String,Vector<String>> expsOfGroup = new LinkedHashMap<String,Vector<String>>();
    if (analysisModule_.getGroupNames()!=null && analysisModule_.getGroupNames().size()>0){
      for (int i=0; i!=analysisModule_.getGroupNames().size(); i++){
        String group = analysisModule_.getGroupNames().get(i);
        expsOfGroup.put(group, analysisModule_.getExpsOfGroup(group));
      }
      //when there are no groups defined -> make one "undefined" group where all experiments belong to
    }else{
      String group = "undefined";
      expsOfGroup.put(group, analysisModule_.getExpNamesInSequence());
    }
    return expsOfGroup;
  }
  
  public ComparativeResultsLookup getComparativeResultsLookup(){
    return analysisModule_;
  }
  
  /**
   * checks whether it is appropriate to initiate an MS1 or MS2 view
   * @param params the data to display
   * @return true when an MS2 view is displayed
   */
  private boolean initMS1OrMS2View(LipidParameterSet params){
    boolean showMS2 = false;
    if (this.displaysMs2_ && params instanceof LipidomicsMSnSet){
      showMS2 = this.showMs2(params,true);
    }else{
      if (!this.displayTolerancePanel_.getLockMzRange().isSelected() || this.viewer_==null || lockRangeUpdateRequired_){
        this.initANewViewer(params);
        lockRangeUpdateRequired_ = false;
      }else {
        updateViewForLockMz(params);
      }        
    }
    return showMS2;
  }
  
  /**
   * this updates the coloring in the 3D viewer, and the 2D chromatogram, when a locked m/z range is used 
   * @param params the selected LipidParameterSet - can be null
   */
  private void updateViewForLockMz(LipidParameterSet params){
    Vector<Vector<CgProbe>> allProbes = getAllProbesFromParams(params,null);
    viewer_.setPaintableProbes(allProbes);
    Vector<CgProbe> storedProbes = allProbes.get(0);
    Vector<CgProbe> selectedProbes = allProbes.get(1);
    this.viewer_.setTheShowSelectedWasOn(false);   
    if (params!=null && this.displayTolerancePanel_.getShow2D().isSelected()) {
      int charge = 1;
      if (params.ProbeCount()>0)
        charge = params.Probe(0).Charge;
      else if (params.getCharge()!=null&&params.getCharge()>1)
        charge = params.getCharge();
      float currentIsotopicMass = params.Mz[0]+(LipidomicsConstants.getNeutronMass()*Integer.parseInt((String)this.isotope_.getSelectedItem())/(float)charge);
      viewer_.setCurrent2DMzRange(currentIsotopicMass-params.LowerMzBand, currentIsotopicMass+params.UpperMzBand);
      float startFloat = currentIsotopicMass-params.LowerMzBand*2;
      float stopFloat = currentIsotopicMass+params.LowerMzBand*2;
      try {
        String[] rawLines = reader_.getRawLines(startFloat, stopFloat, result_.getMsLevels().get(currentSelectedSheet_));
        Hashtable<Integer,Float> rtTimes = reader_.getRetentionTimes();
        Lipidomics2DPainter l2DPainter = new Lipidomics2DPainter(analyzer_,rawLines, rtTimes, reader_.getRetentionTimes(),
            startFloat,stopFloat,reader_.getMultiplicationFactorForInt_()/reader_.getLowestResolution_(),params_.LowerMzBand*2,this,
            MSMapViewer.DISPLAY_TIME_MINUTES,currentIsotopicMass-params.LowerMzBand, currentIsotopicMass+params.UpperMzBand, false,
            storedProbes,selectedProbes,Integer.parseInt((String)this.isotope_.getSelectedItem()),charge,result_.getMsLevels().get(currentSelectedSheet_),
            shotgunIsDisplayed_);
        l2DPainter.preChromatogramExtraxtion(currentIsotopicMass-params.LowerMzBand, currentIsotopicMass+params.UpperMzBand);
        remove2DPainter();
        if (l2dPanelLayout_.getLayoutComponent(BorderLayout.CENTER)!=null && !(l2dPanelLayout_.getLayoutComponent(BorderLayout.CENTER) instanceof Lipidomics2DPainter))
          l2dPanel_.remove(l2dPanelLayout_.getLayoutComponent(BorderLayout.CENTER));
        l2DPainter_ = l2DPainter;
        l2DPainter_.setBackground(Color.WHITE);
        l2dPanel_.add(l2DPainter,BorderLayout.CENTER);
        l2dPanel_.repaint();
        l2dPanel_.invalidate();
        l2dPanel_.updateUI();
//        this.displayPanel_.invalidate();
//        this.majorSplitPane_.repaint();
//        System.out.println("111111111111");

      }
      catch (CgException e) {
        new WarningMessage(new JFrame(), "ERROR", e.getMessage());
        e.printStackTrace();
      }

    }
    this.viewer_.repaintColors();
  }
  
  public boolean isMS2Showing(){
    return this.displaysMs2_;
  }
  
  /**
   * switches between locked and dynamic m/z range
   */
  private void changeLockMzRange() {
    boolean selected = displayTolerancePanel_.getLockMzRange().isSelected();
    displayTolerancePanel_.getDisplayMinusTolerance().setEnabled(!selected);
    displayTolerancePanel_.getDisplayPlusTolerance().setEnabled(!selected);
    displayTolerancePanel_.getDisplayMzStart().setEnabled(selected);
    displayTolerancePanel_.getDisplayMzStop().setEnabled(selected);
    if (displayTolerancePanel_.getLockMzRange().isSelected() && params_!=null && (displayTolerancePanel_.getDisplayMzStart().getText()==null || displayTolerancePanel_.getDisplayMzStart().getText().length()==0)) {
      float currentIsotopicMass = params_.Mz[0];
      displayTolerancePanel_.getDisplayMzStart().setText(Calculator.FormatNumberToString((double)(currentIsotopicMass-Float.parseFloat(this.displayTolerancePanel_.getDisplayMinusTolerance().getText())),2d));
      displayTolerancePanel_.getDisplayMzStop().setText(Calculator.FormatNumberToString((double)(currentIsotopicMass+Float.parseFloat(this.displayTolerancePanel_.getDisplayPlusTolerance().getText())),2d));
    }
    if (displayTolerancePanel_.getLockMzRange().isSelected())
      lockRangeUpdateRequired_ = true;
  }
  
  
  private class HeatMapBuilder implements Runnable
  {
  	private Hashtable<String,ResultDisplaySettings> displaySettingHash_;
    private String molGroup_;
    private Hashtable<String,JPanel> jPanels_;
    
    private HeatMapBuilder(Hashtable<String,ResultDisplaySettings> displaySettingHash, String molGroup, Hashtable<String,JPanel> jPanels)
    {
    	this.displaySettingHash_ = displaySettingHash;
    	this.molGroup_ = molGroup;
    	this.jPanels_ = jPanels;
    }

		@Override
		public void run()
		{
			Hashtable<String,Hashtable<String,Integer>> corrTypeISLookup = analysisModule_.getCorrectionTypeISLookup();
			Hashtable<String,Hashtable<String,Integer>> corrTypeESLookup = analysisModule_.getCorrectionTypeESLookup();
			
			String chosenMolGroup = quantSettingsPanel_.getChosenClassLookup(molGroup_);
			
			JPanel aResultsViewPanel = new JPanel(new BorderLayout());
      JTabbedPane resultsViewTabs= new JTabbedPane();
      aResultsViewPanel.add(resultsViewTabs,BorderLayout.CENTER);
      Hashtable<String,Hashtable<String,ResultCompVO>> resultsOfOneGroup = analysisModule_.getResults().get(molGroup_);
      Vector<String> molNames = analysisModule_.getAllMoleculeNames().get(molGroup_);
      if (molNames.isEmpty()) return; //do not add empty heatmaps
      Hashtable<String,Integer> isLookup = new Hashtable<String,Integer> ();
      Hashtable<String,Integer> esLookup = new Hashtable<String,Integer> ();
      if (corrTypeISLookup.containsKey(chosenMolGroup))
        isLookup = corrTypeISLookup.get(chosenMolGroup);
      if (corrTypeESLookup.containsKey(chosenMolGroup))
        esLookup = corrTypeESLookup.get(chosenMolGroup);
      
//      JPanel aPanel = new JPanel();
      
//      aPanel.setLayout(new BorderLayout());
      boolean hasAbs = jButtonResultAbsQuant_.getText().equalsIgnoreCase("Remove absolute settings");
      ResultDisplaySettings displaySettings = new ResultDisplaySettings(
      		analysisModule_.getISAvailability().get(chosenMolGroup),
      		analysisModule_.getESAvailability().get(chosenMolGroup),isLookup,esLookup,hasAbs,quantSettingsPanel_);
      ResultSelectionSettings selectionSettings = new ResultSelectionSettings("Select molecules to be displayed in the heatmap",molNames,true);
      ResultSelectionSettings combinedChartSettings = new ResultSelectionSettings("Select molecules for a combined bar chart",molNames,false);
      
      HeatMapDrawing drawing = new HeatMapDrawing(resultsOfOneGroup,analysisModule_.getExpNamesInSequence(),molNames, isLookup,esLookup, resultStatus_,LipidDataAnalyzer.this,molGroup_,null,
          displaySettings,selectionSettings,combinedChartSettings,exportSettings_,analysisModule_);
      
      JScrollPane scrollPane = new JScrollPane(drawing);
      scrollPane.getViewport().getView().setBackground(Color.GRAY);
      int[] widthAndHeight = getScrollPaneWidthAndHeight(drawing);
      scrollPane.setPreferredSize(new Dimension(widthAndHeight[0], widthAndHeight[1]));
      heatmaps_.put(molGroup_, drawing);
      
      JSplitPane splitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT, scrollPane, new JScrollPane(drawing.getSettingsPanel()));
			splitPane.setOneTouchExpandable(true);
			splitPane.setDividerLocation(600);
      
      resultsViewTabs.addTab("Heatmap", splitPane);
      resultsViewTabs.setToolTipTextAt(0, TooltipTexts.TABS_RESULTS_HEATMAP+molGroup_+"</html>");
      JPanel barChartPanel = new JPanel();
      resultsViewTabs.addTab("Bar-chart", barChartPanel);
      resultsViewTabs.setToolTipTextAt(1, TooltipTexts.TABS_RESULTS_BARCHART+molGroup_+"</html>");
      if (groupsPanel_.getGroups().size()>0){
        Hashtable<String,Hashtable<String,ResultCompVO>> groupedResultsOfOneGroup = analysisModule_.getGroupedResults().get(molGroup_);
        
        HeatMapDrawing groupDrawing = new HeatMapDrawing(groupedResultsOfOneGroup, groupsPanel_.getGroups(),molNames, isLookup,esLookup, resultStatus_,LipidDataAnalyzer.this,molGroup_,drawing, 
        		displaySettings,selectionSettings,combinedChartSettings,exportSettingsGroup_,analysisModule_);
        
        JScrollPane groupScrollPane = new JScrollPane(groupDrawing);
        groupScrollPane.getViewport().getView().setBackground(Color.GRAY);
        int[] groupWidthAndHeight = getScrollPaneWidthAndHeight(groupDrawing);
        groupScrollPane.setPreferredSize(new Dimension(groupWidthAndHeight[0], groupWidthAndHeight[1]));
        groupHeatmaps_.put(molGroup_, groupDrawing);
        
        JSplitPane groupSplitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT, groupScrollPane, new JScrollPane(groupDrawing.getSettingsPanel()));
        groupSplitPane.setOneTouchExpandable(true);
        groupSplitPane.setDividerLocation(640);
        
        resultsViewTabs.addTab("Group-Heatmap", groupSplitPane);
        resultsViewTabs.setToolTipTextAt(2, TooltipTexts.TABS_RESULTS_HEATMAP_GROUP+molGroup_+"</html>");
        JPanel groupBarChartPanel = new JPanel();
        resultsViewTabs.addTab("Group bar-chart", groupBarChartPanel);
        resultsViewTabs.setToolTipTextAt(3, TooltipTexts.TABS_RESULTS_BARCHART_GROUP+molGroup_+"</html>");
      }
      displaySettingHash_.put(molGroup_, displaySettings);
      molBarCharts_.put(molGroup_, resultsViewTabs);
      jPanels_.put(molGroup_, aResultsViewPanel);
//      resultTabs_.addTab(molGroup_,aResultsViewPanel);
//      resultTabs_.setToolTipTextAt(resultTabs_.indexOfTab(molGroup_), TooltipTexts.TABS_RESULTS_GROUP+molGroup_+"</html>");
		}
  }
  
  
}
  
 