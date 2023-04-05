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
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
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
import java.util.Timer;
import java.util.TimerTask;
import java.util.Vector;

import javax.imageio.ImageIO;
import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.Icon;
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
import javax.swing.JProgressBar;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.ListSelectionModel;
import javax.swing.UIManager;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.table.TableColumnModel;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.TransformerConfigurationException;
import javax.xml.transform.TransformerException;

import uk.ac.ebi.pride.jmztab2.utils.errors.MZTabErrorType.Level;

import org.apache.batik.dom.GenericDOMImplementation;
import org.apache.batik.svggen.SVGGraphics2D;
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
import at.tugraz.genome.lda.exception.RdbWriterException;
import at.tugraz.genome.lda.exception.RulesException;
import at.tugraz.genome.lda.exception.SettingsException;
import at.tugraz.genome.lda.export.LDAExporter;
import at.tugraz.genome.lda.export.OmegaMasslistExporter;
import at.tugraz.genome.lda.export.QuantificationResultExporter;
import at.tugraz.genome.lda.export.vos.AnalyteOmegaInfoVO;
import at.tugraz.genome.lda.interfaces.ColorChangeListener;
import at.tugraz.genome.lda.listeners.AnnotationThresholdListener;
import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
import at.tugraz.genome.lda.msn.MSnAnalyzer;
import at.tugraz.genome.lda.msn.hydroxy.parser.HydroxyEncoding;
import at.tugraz.genome.lda.mztab.MztabUtils;
import at.tugraz.genome.lda.mztab.SmallMztabMolecule;
import at.tugraz.genome.lda.parser.LDAResultReader;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.quantification.LipidomicsAnalyzer;
import at.tugraz.genome.lda.quantification.QuantificationResult;
import at.tugraz.genome.lda.swing.AbsoluteQuantSettingsPanel;
import at.tugraz.genome.lda.swing.BarChartPainter;
import at.tugraz.genome.lda.swing.BatchQuantificationTable;
import at.tugraz.genome.lda.swing.BatchQuantificationTableModel;
import at.tugraz.genome.lda.swing.ClassesOverviewPanel;
import at.tugraz.genome.lda.swing.ColorChooserDialog;
import at.tugraz.genome.lda.swing.CutoffSettingsPanel;
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
import at.tugraz.genome.lda.swing.OmegaExportDialog;
import at.tugraz.genome.lda.swing.RangeColor;
import at.tugraz.genome.lda.swing.RecalculateMSnDialog;
import at.tugraz.genome.lda.swing.ResultDisplaySettings;
import at.tugraz.genome.lda.swing.ResultSelectionSettings;
import at.tugraz.genome.lda.swing.RuleDefinitionInterface;
import at.tugraz.genome.lda.swing.SpectrumUpdateListener;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.verifier.DoubleVerifier;
import at.tugraz.genome.lda.verifier.IntegerMaxVerifier;
import at.tugraz.genome.lda.vos.AbsoluteSettingsVO;
import at.tugraz.genome.lda.vos.AddAnalyteVO;
import at.tugraz.genome.lda.vos.AutoAnalyteAddVO;
import at.tugraz.genome.lda.vos.IntegerStringVO;
import at.tugraz.genome.lda.vos.IsotopicLabelVO;
import at.tugraz.genome.lda.vos.QuantVO;
import at.tugraz.genome.lda.vos.RawQuantificationPairVO;
import at.tugraz.genome.lda.vos.ResultAreaVO;
import at.tugraz.genome.lda.vos.ResultCompVO;
import at.tugraz.genome.lda.vos.ResultDisplaySettingsVO;
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

import com.sun.j3d.utils.applet.MainFrame;

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
  
  private static MainFrame frame_;
  private JFileChooser mzXMLFileChooser;
  private JFileChooser mzXMLDirChooser;
  private JTabbedPane mainTabs;
  private JPanel singleQuantMenu_;
  private JPanel batchQuantMenu_;
  private JPanel resultsMenu_;
  private JPanel resultsPanel_;
  private JPanel settingsPanel_;
  private JPanel helpPanel_;
  private JPanel aboutPanel_;
  private JPanel displayTopMenu;
  private JLabel mzXMLLabel;
  private JTextField selectedMzxmlFile;
  private JTextField selectedMzxmlDirectory_;
  private JButton jButtonMxXMLOpen;
  private JButton jButtonMxXMLDirOpen_;
  private QuantificationThread quantThread_;
  private BatchQuantThread batchQuantThread_;
  private RawToMzxmlThread rawmzThread_;
  private MzxmlToChromThread mzToChromThread_;
  private Timer timer_;
  private JLabel quantifyingLabel_;
  private JLabel quantifyingBatchLabel_;
  
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
  private JScrollPane analysisTablePane;
  private JPanel analysisSelectionTablePanel_;
  private JTable resultFilesDisplayTable;
  private ListSelectionModel resultListSelectionModel;
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

  
  private JLabel chromLabel;
  private JTextField selectedChromFile;
  private JButton jButtonChromOpen;
  private JTextField selectedQuantFile;
  private JTextField selectedQuantDir_;
  private JButton jButtonQuantOpen;
  private JButton jBatchQuantOpen_;
  private JButton startQuantification;
  private JButton startBatchQuantification_;
  private JFileChooser chromFileChooser_;
  private JFileChooser quantFileChooser_;
  private JFileChooser quantDirChooser_;
  private JTextField singleTimeMinusTol_;
  private JTextField singleTimePlusTol_;
  private JTextField singleCutoff_;
  private JTextField singleRTShift_;
  private JCheckBox isoValidation_;
  private JTextField amountOfIsotopes_;
  private JCheckBox isoBatchValidation_;
  private JTextField amountOfBatchIsotopes_;
  private JCheckBox searchUnknownTime_;
  private JCheckBox searchUnknownBatchTime_;
  private JTextField amountOfMatchingSearchIsotopes_;
  private JTextField amountOfMatchingBatchSearchIsotopes_;
  private JTextField nrProcessors_;
  private JTextField nrProcessorsBatch_;
  /** the ion mode for the Alex123 searches*/
  private JComboBox<String> ionMode_;
  /** the ion mode for the Alex123 batch searches*/
  private JComboBox<String> ionModeBatch_;
  /** the ion mode to set Alex123 order in display results*/
  private JComboBox<String> ionModeOrder_;

  private JTextField batchTimeMinusTol_;
  private JTextField batchTimePlusTol_;
  private JTextField batchCutoff_;
  private JTextField batchRTShift_;

  
  private JTextField selectedResultFile;
  private JButton jButtonResultOpen;
  private JButton startDisplay;
  private JFileChooser resultFileChooser_;
  private LipidomicsJTable displayTable;
  private BatchQuantificationTable batchQuantTable_;
  private BatchQuantificationTableModel batchQuantTableModel_;
  private ListSelectionModel listSelectionModel;
  private JPanel selectionPane;
  private JPanel tableContainer;
  private JPanel tablePanel_;
  private JComboBox<String> selectedSheet_;
  private JPanel quantifyingPanel_;
  private JPanel quantifyingBatchPanel_;  
  private JTextField displayMinusTolerance_;
  private JTextField displayPlusTolerance_;
  /** text field to set a permanent m/z start value for the 3D viewer*/
  private JTextField displayMzStart_;
  /** text field to set a permanent m/z stop value for the 3D viewer*/
  private JTextField displayMzStop_;
  private JTextField displayRtStart_;
  private JTextField displayRtStop_;

  private JCheckBox show2D_;
  /** show the names in MSn style in the display results or not*/
  private JCheckBox showMSnNames_;
  
  /** this check box decides whether a static (locked) range is used, or it is adapted to the selected analyte*/
  private JCheckBox lockMzRange_;

  private Lipidomics2DPainter l2DPainter_;
  private Lipidomics2DPainter spectrumPainter_;
  private LipidomicsItemListener changeIsotopeListener_;
  
  private JPanel displayPanel_;
  private JSplitPane majorSplitPane_;
  private JSplitPane topSplitPane_;
  private RuleDefinitionInterface userInterface;
  private JScrollPane tablePane;
  private QuantificationResult result_;
  private Hashtable<Integer,Integer> resultPositionToOriginalLoopkup_ = new Hashtable<Integer,Integer>();
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
  private JProgressBar progressBar_;
  private JProgressBar progressBatchBar_;
  private JLabel spinnerLabel_;
  private JLabel spinnerBatchLabel_;
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
  
  /** action for setting a locked m/z range for the 3D viewer*/
  private final static String CHANGE_LOCK_MZ_RANGE = "lockMzRange";
  /** the displayed label next to check box for activating the lock mass range*/
  private final static String DISPLAY_LOCK_MZ_TEXT = "Lock m/z range ";
  
  /** Global finalButtonSection to delete it in the makeDisplayRemoveOperations */
  private JPanel finalButtonSection_;
  
  /** Object for the Rule Definition Interface */
  private RuleDefinitionInterface msnUserInterfaceObject_;
  
  private final static String DEFAULT_ANNOTATION_CUTOFF = "5";
  
  private final static Font SELECT_FIELD_FONT = new Font("Helvetica",Font.PLAIN,10);
  private final static Font SMALL_FONT = new Font("Arial",Font.PLAIN,9);
  
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
  /** the export dialog for omega retention time map*/
//  private OmegaExportDialog  omegaExport_ = null;

  
  public LipidDataAnalyzer(){
    this.createDisplayTopMenu();
    this.createSingleQuantMenu();
    this.createBatchQuantMenu();
    this.createResultsMenu();
    this.initL2dPanel();
    displaysMs2_ = false;
    shotgunIsDisplayed_ = false;

    JPanel displayTolerancePanel = new JPanel();
    displayTolerancePanel.setLayout(new GridBagLayout());
    JLabel diplayTolMinus = new JLabel("- m/z: ");
    diplayTolMinus.setFont(SMALL_FONT);
    diplayTolMinus.setToolTipText(TooltipTexts.DISPLAY_MZ_MINUS);
    displayTolerancePanel.add(diplayTolMinus,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 1, 0, 0), 0, 0));
    JLabel diplayTolPlus = new JLabel("+ m/z: ");
    diplayTolPlus.setFont(SMALL_FONT);
    diplayTolPlus.setToolTipText(TooltipTexts.DISPLAY_MZ_PLUS);
    displayTolerancePanel.add(diplayTolPlus,new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 1, 0, 0), 0, 0));
    displayMinusTolerance_ = new JTextField(2);
    displayMinusTolerance_.setFont(SMALL_FONT);
    displayMinusTolerance_.setText("1.5");
    displayMinusTolerance_.setHorizontalAlignment(JTextField.RIGHT);
    displayMinusTolerance_.setToolTipText(TooltipTexts.DISPLAY_MZ_MINUS);
    displayMinusTolerance_.setInputVerifier(new DoubleVerifier());
    displayTolerancePanel.add(displayMinusTolerance_,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 1, 0, 0), 0, 0));    
    displayPlusTolerance_ = new JTextField(2);
    displayPlusTolerance_.setFont(SMALL_FONT);
    displayPlusTolerance_.setText("2.5");
    displayPlusTolerance_.setHorizontalAlignment(JTextField.RIGHT);
    displayPlusTolerance_.setToolTipText(TooltipTexts.DISPLAY_MZ_PLUS);
    displayPlusTolerance_.setInputVerifier(new DoubleVerifier());
    displayTolerancePanel.add(displayPlusTolerance_,new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 1, 0, 0), 0, 0));
    JLabel diplayTolUnit1 = new JLabel("[Da]");
    diplayTolUnit1.setFont(SMALL_FONT);
    displayTolerancePanel.add(diplayTolUnit1,new GridBagConstraints(4, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 1, 0, 0), 0, 0));
    JLabel rtTolMinusText = new JLabel("Start: ");
    rtTolMinusText.setFont(SMALL_FONT);
    rtTolMinusText.setToolTipText(TooltipTexts.DISPLAY_RT_START);
    displayTolerancePanel.add(rtTolMinusText,new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 1, 0, 0), 0, 0));
    JLabel rtTolPlusText = new JLabel("Stop: ");
    rtTolPlusText.setFont(SMALL_FONT);
    rtTolPlusText.setToolTipText(TooltipTexts.DISPLAY_RT_STOP);
    displayTolerancePanel.add(rtTolPlusText,new GridBagConstraints(2, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 1, 0, 0), 0, 0));
    displayRtStart_ = new JTextField(2);
    displayRtStart_.setFont(SMALL_FONT);
    displayRtStart_.setText("");
    displayRtStart_.setHorizontalAlignment(JTextField.RIGHT);
    displayRtStart_.setToolTipText(TooltipTexts.DISPLAY_RT_START);
    displayRtStart_.setInputVerifier(new DoubleVerifier());
    displayTolerancePanel.add(displayRtStart_,new GridBagConstraints(1, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 1, 0, 0), 0, 0));
    displayRtStop_ = new JTextField(2);
    displayRtStop_.setFont(SMALL_FONT);
    displayRtStop_.setText("");
    displayRtStop_.setHorizontalAlignment(JTextField.RIGHT);
    displayRtStop_.setToolTipText(TooltipTexts.DISPLAY_RT_STOP);
    displayRtStop_.setInputVerifier(new DoubleVerifier());
    displayTolerancePanel.add(displayRtStop_,new GridBagConstraints(3, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 1, 0, 0), 0, 0));
    JLabel diplayTolUnit2 = new JLabel("[min]");
    diplayTolUnit2.setFont(SMALL_FONT);
//    diplayTolUnit1.setToolTipText(TooltipTexts.DISPLAY_MZ_MINUS);
    displayTolerancePanel.add(diplayTolUnit2,new GridBagConstraints(4, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 1, 0, 0), 0, 0));
    
    JPanel lockMxSelectionPanel = new JPanel();
    lockMzRange_ = new JCheckBox();
    lockMzRange_.setFont(SMALL_FONT);
    lockMzRange_.addItemListener(new LipidomicsItemListener(CHANGE_LOCK_MZ_RANGE));
    lockMzRange_.setToolTipText(TooltipTexts.DISPLAY_LOCK_MZ);
    lockMxSelectionPanel.add(lockMzRange_);
    JLabel lockMz = new JLabel(DISPLAY_LOCK_MZ_TEXT);
    lockMz.setToolTipText(TooltipTexts.DISPLAY_LOCK_MZ);
    lockMz.setFont(SMALL_FONT);
    lockMxSelectionPanel.add(lockMz);
    displayMzStart_ = new JTextField(4);
    displayMzStart_.setFont(SMALL_FONT);
    displayMzStart_.setText("");
    displayMzStart_.setHorizontalAlignment(JTextField.RIGHT);
    displayMzStart_.setToolTipText(TooltipTexts.DISPLAY_MZ_MINUS);
    displayMzStart_.setInputVerifier(new DoubleVerifier());
    lockMxSelectionPanel.add(displayMzStart_);
    displayMzStart_.setEnabled(false);
    JLabel mzRangeSign = new JLabel("-");
    lockMxSelectionPanel.add(mzRangeSign);
    displayMzStop_ = new JTextField(4);
    displayMzStop_.setFont(SMALL_FONT);
    displayMzStop_.setText("");
    displayMzStop_.setHorizontalAlignment(JTextField.RIGHT);
    displayMzStop_.setToolTipText(TooltipTexts.DISPLAY_MZ_MINUS);
    displayMzStop_.setInputVerifier(new DoubleVerifier());
    displayMzStop_.setEnabled(false);
    lockMxSelectionPanel.add(displayMzStop_);
    
    displayTolerancePanel.add(lockMxSelectionPanel ,new GridBagConstraints(0, 2, 6, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 1, 0, 0), 0, 0));

    
    JButton quantTolUpdate = new JButton("Update");
    quantTolUpdate.addActionListener(this);
    quantTolUpdate.setFont(quantTolUpdate.getFont().deriveFont(10f));
    quantTolUpdate.setMargin(new Insets(1,5,1,5));
    quantTolUpdate.setActionCommand("updateQuantTolOfCurrentlySelected");
    quantTolUpdate.setToolTipText(TooltipTexts.DISPLAY_UPDATE);
    displayTolerancePanel.add(quantTolUpdate,new GridBagConstraints(5, 0, 1, 2, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    JPanel showOptionsPanel = new JPanel();
    displayTolerancePanel.add(showOptionsPanel,new GridBagConstraints(0, 3, 6, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    
    showMSnNames_ = new JCheckBox();
    showMSnNames_.addItemListener(new LipidomicsItemListener("showMSnNames"));
    showMSnNames_.setToolTipText(TooltipTexts.DISPLAY_SHOW_MSN);
    showOptionsPanel.add(showMSnNames_);
    JLabel showMSn = new JLabel("Show MSn");
    showMSn.setToolTipText(TooltipTexts.DISPLAY_SHOW_MSN);
    showOptionsPanel.add(showMSn);
    
    show2D_ = new JCheckBox();
    show2D_.setSelected(true);
    show2D_.addItemListener(new LipidomicsItemListener("show2dChanged"));
    show2D_.setToolTipText(TooltipTexts.DISPLAY_SHOW_2D);
    showOptionsPanel.add(show2D_);
    JLabel show2d = new JLabel("Show 2D-View");
    show2d.setToolTipText(TooltipTexts.DISPLAY_SHOW_2D);
    showOptionsPanel.add(show2d);

    
    selectionPane = new JPanel();
    selectionPane.setLayout(new BoxLayout(selectionPane, BoxLayout.LINE_AXIS));
    Vector<LipidParameterSet> dummy = new Vector<LipidParameterSet>();
    

    displayTable = new LipidomicsJTable(new LipidomicsTableModel(dummy,dummy,false,false),new LipidomicsTableCellRenderer(),false,0,true,this);
    listSelectionModel = displayTable.getSelectionModel();
    displayTable.setSelectionModel(listSelectionModel);
    listSelectionModel.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
    
    JPanel listContainer = new JPanel(new GridLayout(1,1));
    tablePane = new JScrollPane(displayTable);

    
    tablePanel_ = new JPanel();
    tablePanel_.setLayout(new BorderLayout());
    selectedSheet_ = new JComboBox<String>();
    tablePanel_.add(selectedSheet_,BorderLayout.NORTH);
    tablePanel_.add(tablePane,BorderLayout.CENTER);
    tablePanel_.add(displayTolerancePanel,BorderLayout.SOUTH);

    
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

    JPanel singleQuantificationPanel_ = new JPanel();
    singleQuantificationPanel_.setLayout(new BorderLayout());
    JPanel batchQuantificationPanel_ = new JPanel();
    batchQuantificationPanel_.setLayout(new BorderLayout());
    resultsPanel_ = new JPanel();
    resultsPanel_.setLayout(new BorderLayout());
    settingsPanel_ = new JPanel();
    initSettingsPanel();
    helpPanel_ = new JPanel();
    initHelpPanel();
    aboutPanel_ = new JPanel();
    initAboutPanel();

    
    mainTabs.addTab("Quantitation", singleQuantificationPanel_);
    mainTabs.setToolTipTextAt(0, TooltipTexts.TABS_MAIN_QUANTITATION);
    mainTabs.addTab("Batch Quantitation", batchQuantificationPanel_);
    mainTabs.setToolTipTextAt(1, TooltipTexts.TABS_MAIN_BATCH);
    mainTabs.addTab("Statistical Analysis", resultsPanel_);
    mainTabs.setToolTipTextAt(2, TooltipTexts.TABS_MAIN_STATISTICS);
    mainTabs.addTab("Display Results", displayPanel_);
    mainTabs.setToolTipTextAt(3, TooltipTexts.TABS_MAIN_DISPLAY);
    mainTabs.addTab("Settings", settingsPanel_);
    mainTabs.setToolTipTextAt(4, TooltipTexts.TABS_MAIN_SETTINGS);
    mainTabs.addTab("Help", helpPanel_);
    mainTabs.setToolTipTextAt(5, TooltipTexts.TABS_MAIN_HELP);
    mainTabs.addTab("About", aboutPanel_);
    mainTabs.setToolTipTextAt(6, TooltipTexts.TABS_MAIN_ABOUT);
    LicenseChangeListener licenseListener = new LicenseChangeListener();
    //mainTabs.addMouseListener(licenseListener);
    mainTabs.addChangeListener(licenseListener);
    mainTabs.setSelectedIndex(0);
    singleQuantificationPanel_.add(singleQuantMenu_);
    batchQuantificationPanel_.add(batchQuantMenu_);
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
    linkAddress = "http://genome.tugraz.at/lda2/data/review/DataDescription.pdf";
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
    text = new JLabel("J\u00fcrgen Hartler, Alexander Triebl, Martin Tr\u00f6tzm\u00fcller, Andreas Ziegl");
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
    text = new JLabel("Copyright \u00A9 2023 J\u00fcrgen Hartler, Andreas Ziegl, Gerhard G Thallinger, Leonida M Lamp");
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
    selectedChromFile.setToolTipText(TooltipTexts.DISPLAY_OPEN_CHROM);
    displayTopMenu.add(selectedChromFile,new GridBagConstraints(0, 0, 6, 1, 0.0, 0.0
      ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    jButtonChromOpen = new JButton("Open Chrom");
    jButtonChromOpen.addActionListener(this);
    jButtonChromOpen.setActionCommand("showChromFileChooser");
    jButtonChromOpen.setToolTipText(TooltipTexts.DISPLAY_OPEN_CHROM);
    displayTopMenu.add(jButtonChromOpen,new GridBagConstraints(7, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    this.selectedResultFile = new JTextField(62);
    selectedResultFile.setToolTipText(TooltipTexts.DISPLAY_OPEN_RESULT);
    displayTopMenu.add(selectedResultFile,new GridBagConstraints(0, 1, 6, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    this.jButtonResultOpen = new JButton("Open Result");
    jButtonResultOpen.addActionListener(this);
    jButtonResultOpen.setActionCommand("showResultChooser");
    jButtonResultOpen.setToolTipText(TooltipTexts.DISPLAY_OPEN_RESULT);
    displayTopMenu.add(jButtonResultOpen,new GridBagConstraints(7, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    startDisplay = new JButton("Start Display");
    startDisplay.addActionListener(this);
    startDisplay.setActionCommand("startDisplay");
    startDisplay.setToolTipText(TooltipTexts.DISPLAY_START);
    displayTopMenu.add(startDisplay,new GridBagConstraints(8, 0, 1, 2, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0)); 
    
    jButtonResultOpen.setPreferredSize(jButtonChromOpen.getPreferredSize());
  }
  
  private void createBatchQuantMenu(){
    this.batchQuantMenu_ = new JPanel();
    this.batchQuantMenu_.setLayout(new GridBagLayout());
    JPanel selectionPanel = new JPanel();
    selectionPanel.setLayout(new GridBagLayout());
    batchQuantMenu_.add(selectionPanel,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    JLabel batchMzXMLLabel = new JLabel("Raw files: ");
    batchMzXMLLabel.setToolTipText(TooltipTexts.QUANTITATION_BATCH_RAW_FILE);
    selectionPanel.add(batchMzXMLLabel,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    this.selectedMzxmlDirectory_ = new JTextField(62);
    selectedMzxmlDirectory_.setToolTipText(TooltipTexts.QUANTITATION_BATCH_RAW_FILE);
    selectionPanel.add(selectedMzxmlDirectory_,new GridBagConstraints(1, 0, 6, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    jButtonMxXMLDirOpen_ = new JButton("Select");
    jButtonMxXMLDirOpen_.addActionListener(this);
    jButtonMxXMLDirOpen_.setActionCommand("showMzxmlDirChooser");
    jButtonMxXMLDirOpen_.setToolTipText(TooltipTexts.QUANTITATION_BATCH_RAW_FILE);
    selectionPanel.add(jButtonMxXMLDirOpen_,new GridBagConstraints(7, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    JLabel quantChromLabel = new JLabel("Quant. files: ");
    quantChromLabel.setToolTipText(TooltipTexts.QUANTITATION_BATCH_MASS_LIST);
    selectionPanel.add(quantChromLabel,new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    jBatchQuantOpen_ = new JButton("Select");
    jBatchQuantOpen_.addActionListener(this);
    jBatchQuantOpen_.setActionCommand("showQuantDirChooser");
    jBatchQuantOpen_.setToolTipText(TooltipTexts.QUANTITATION_BATCH_MASS_LIST);
    selectionPanel.add(jBatchQuantOpen_,new GridBagConstraints(7, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    this.selectedQuantDir_ = new JTextField(62);
    selectedQuantDir_.setToolTipText(TooltipTexts.QUANTITATION_BATCH_MASS_LIST);
    selectionPanel.add(selectedQuantDir_,new GridBagConstraints(1, 1, 6, 1, 0.0, 0.0
      ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    
    JPanel settingsPanel = new JPanel ();  
    settingsPanel.setLayout(new GridBagLayout());
    batchQuantMenu_.add(settingsPanel,new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(6, 6, 0, 0), 0, 0));

//    JLabel mzTolLabel = new JLabel("m/z-Tolerance: ");
//    settingsPanel.add(mzTolLabel,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
//        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
//    batchMzTol_ = new JTextField(4);
//    batchMzTol_.setText("0.02");
//    batchMzTol_.setHorizontalAlignment(JTextField.RIGHT);
//    settingsPanel.add(batchMzTol_,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
//        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    JLabel timeMinusTolLabel = new JLabel("Time before tol.: ");
    timeMinusTolLabel.setToolTipText(TooltipTexts.QUANTITATION_RET_BEFORE);
    settingsPanel.add(timeMinusTolLabel,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    batchTimeMinusTol_ = new JTextField(4);
    batchTimeMinusTol_.setText("5");
    batchTimeMinusTol_.setHorizontalAlignment(JTextField.RIGHT);
    batchTimeMinusTol_.setToolTipText(TooltipTexts.QUANTITATION_RET_BEFORE);
    settingsPanel.add(batchTimeMinusTol_,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    JLabel timeMinusUnit = new JLabel("min");
    timeMinusUnit.setToolTipText(TooltipTexts.QUANTITATION_RET_BEFORE);
    settingsPanel.add(timeMinusUnit,new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    JLabel timePlusTolLabel = new JLabel("Time after tol.: ");
    timePlusTolLabel.setToolTipText(TooltipTexts.QUANTITATION_RET_AFTER);
    settingsPanel.add(timePlusTolLabel,new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    batchTimePlusTol_ = new JTextField(4);
    batchTimePlusTol_.setText("5");
    batchTimePlusTol_.setHorizontalAlignment(JTextField.RIGHT);
    batchTimePlusTol_.setToolTipText(TooltipTexts.QUANTITATION_RET_AFTER);
    settingsPanel.add(batchTimePlusTol_,new GridBagConstraints(4, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    JLabel timePlusUnit = new JLabel("min");
    timePlusUnit.setToolTipText(TooltipTexts.QUANTITATION_RET_AFTER);
    settingsPanel.add(timePlusUnit,new GridBagConstraints(5, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));

    JLabel cutoffLabel = new JLabel("Rel. base-peak cutoff: ");
    cutoffLabel.setToolTipText(TooltipTexts.QUANTITATION_CUTOFF);
    settingsPanel.add(cutoffLabel,new GridBagConstraints(6, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    batchCutoff_ = new JTextField(4);
    batchCutoff_.setText(LipidomicsConstants.getBasePeakDefaultCutoff());
    batchCutoff_.setHorizontalAlignment(JTextField.RIGHT);
    batchCutoff_.setToolTipText(TooltipTexts.QUANTITATION_CUTOFF);
    settingsPanel.add(batchCutoff_,new GridBagConstraints(7, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    JLabel cutoffUnit = new JLabel("\u2030");
    cutoffUnit.setToolTipText(TooltipTexts.QUANTITATION_CUTOFF);
    settingsPanel.add(cutoffUnit,new GridBagConstraints(8, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    JLabel rtShiftLabel = new JLabel("RT-shift: ");
    rtShiftLabel.setToolTipText(TooltipTexts.QUANTIFICATION_RET_SHIFT);
    settingsPanel.add(rtShiftLabel,new GridBagConstraints(9, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    batchRTShift_ = new JTextField(4);
    batchRTShift_.setText("0.0");
    batchRTShift_.setHorizontalAlignment(JTextField.RIGHT);
    batchRTShift_.setToolTipText(TooltipTexts.QUANTIFICATION_RET_SHIFT);
    settingsPanel.add(batchRTShift_,new GridBagConstraints(10, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    JLabel rtShiftUnit = new JLabel("min");
    rtShiftUnit.setToolTipText(TooltipTexts.QUANTIFICATION_RET_SHIFT);
    settingsPanel.add(rtShiftUnit,new GridBagConstraints(11, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));

    
    JPanel settingsPanel2 = new JPanel ();
    settingsPanel2.setLayout(new GridBagLayout());
    batchQuantMenu_.add(settingsPanel2,new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(6, 6, 0, 0), 0, 0));
    isoBatchValidation_ = new JCheckBox();
    isoBatchValidation_.setSelected(true);
    isoBatchValidation_.addItemListener(new LipidomicsItemListener("endisableAmountOfBatchIsotopes"));
    isoBatchValidation_.setToolTipText(TooltipTexts.QUANTITATION_ISOTOPES);
    settingsPanel2.add(isoBatchValidation_,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    JLabel allowIsotopicValidation = new JLabel("Isotopic quantitation of ");
    allowIsotopicValidation.setToolTipText(TooltipTexts.QUANTITATION_ISOTOPES);
    settingsPanel2.add(allowIsotopicValidation,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    amountOfBatchIsotopes_ = new JTextField(4);
    amountOfBatchIsotopes_.setText("2");
    amountOfBatchIsotopes_.setHorizontalAlignment(JTextField.RIGHT);
    amountOfBatchIsotopes_.setToolTipText(TooltipTexts.QUANTITATION_ISOTOPES);
    settingsPanel2.add(amountOfBatchIsotopes_,new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    JLabel isotopePeaks = new JLabel("isotopes where ");
    isotopePeaks.setToolTipText(TooltipTexts.QUANTITATION_ISOTOPES);
    settingsPanel2.add(isotopePeaks,new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    amountOfMatchingBatchSearchIsotopes_ = new JTextField(4);
    amountOfMatchingBatchSearchIsotopes_.setText("1");
    amountOfMatchingBatchSearchIsotopes_.setHorizontalAlignment(JTextField.RIGHT);
//    amountOfMatchingBatchSearchIsotopes_.setEnabled(false);
    amountOfMatchingBatchSearchIsotopes_.setToolTipText(TooltipTexts.QUANTITATION_ISOTOPES);
    settingsPanel2.add(amountOfMatchingBatchSearchIsotopes_,new GridBagConstraints(4, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    JLabel isotopeMatchPeaks = new JLabel("isotopic peak(s) have to match");
    isotopeMatchPeaks.setToolTipText(TooltipTexts.QUANTITATION_ISOTOPES);
    settingsPanel2.add(isotopeMatchPeaks,new GridBagConstraints(5, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));


    JPanel settingsPanel3 = new JPanel ();
    settingsPanel3.setLayout(new GridBagLayout());
    batchQuantMenu_.add(settingsPanel3,new GridBagConstraints(0, 3, 1, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(6, 6, 0, 0), 0, 0));
    searchUnknownBatchTime_ = new JCheckBox();
    searchUnknownBatchTime_.setSelected(true);
    searchUnknownBatchTime_.setToolTipText(TooltipTexts.QUANTITATION_RET_UNKNOWN);
    settingsPanel3.add(searchUnknownBatchTime_,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    JLabel findUnknownMolecules = new JLabel("Find molecules where retention time is unknown");
    findUnknownMolecules.setToolTipText(TooltipTexts.QUANTITATION_RET_UNKNOWN);
    settingsPanel3.add(findUnknownMolecules,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));

    JLabel maxProcessors = new JLabel("Processors for quantitation:");
    maxProcessors.setToolTipText(TooltipTexts.QUANTITATION_PROCESSORS);
    settingsPanel3.add(maxProcessors,new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    nrProcessorsBatch_ = new JTextField(3);
    nrProcessorsBatch_.setText(String.valueOf(this.getAmountOfProcessorsPreferred()));
    nrProcessorsBatch_.setHorizontalAlignment(JTextField.RIGHT);
    nrProcessorsBatch_.setToolTipText(TooltipTexts.QUANTITATION_PROCESSORS);
    nrProcessorsBatch_.setInputVerifier(new IntegerMaxVerifier(true,1,getMaxProcessors()));
    settingsPanel3.add(nrProcessorsBatch_,new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));

    int yPositionCounterSettingsPanels = 4;
    if (Settings.useAlex()){
      JPanel settingsPanel4 = new JPanel();
      settingsPanel4.setLayout(new GridBagLayout());
      batchQuantMenu_.add(settingsPanel4,new GridBagConstraints(0, yPositionCounterSettingsPanels, 1, 1, 0.0, 0.0
          ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(6, 6, 0, 0), 0, 0));
      JLabel ionModeLabel = new JLabel("Ion mode:");
      ionModeLabel.setToolTipText(TooltipTexts.QUANTITATION_ION_MODE);
      settingsPanel4.add(ionModeLabel,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
          ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
      ionModeBatch_ = new JComboBox<String>();
      ionModeBatch_.addItem("+");
      ionModeBatch_.addItem("-");
      ionModeBatch_.setToolTipText(TooltipTexts.QUANTITATION_ION_MODE);
      ionModeBatch_.setFont(SELECT_FIELD_FONT);
      settingsPanel4.add(ionModeBatch_,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
          ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
      JLabel ionModeInfo = new JLabel("(selection required for Alex123 target lists)");
      ionModeInfo.setToolTipText(TooltipTexts.QUANTITATION_ION_MODE);
      settingsPanel4.add(ionModeInfo,new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
      yPositionCounterSettingsPanels++;
    }
    
    startBatchQuantification_ = new JButton("Start Quantitation");
    startBatchQuantification_.addActionListener(this);
    startBatchQuantification_.setActionCommand("startBatchQuantification");
    startBatchQuantification_.setToolTipText(TooltipTexts.QUANTITATION_BATCH_START);
    batchQuantMenu_.add(startBatchQuantification_,new GridBagConstraints(0, yPositionCounterSettingsPanels, 1, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(6, 6, 0, 0), 0, 0));
    yPositionCounterSettingsPanels++;
    quantifyingBatchPanel_ = new JPanel ();  
    quantifyingBatchPanel_.setLayout(new BorderLayout());
    batchQuantMenu_.add(quantifyingBatchPanel_,new GridBagConstraints(0, yPositionCounterSettingsPanels, 1, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(6, 6, 0, 0), 0, 0));
    yPositionCounterSettingsPanels++;
    
    batchQuantTableModel_ = new BatchQuantificationTableModel();
    batchQuantTable_ = new BatchQuantificationTable(batchQuantTableModel_);

    TableColumnModel tcm = batchQuantTable_.getColumnModel();
    tcm.getColumn(0).setPreferredWidth(15);
    tcm.getColumn(0).setMaxWidth(15);
    tcm.getColumn(1).setPreferredWidth(170);
    tcm.getColumn(2).setPreferredWidth(200);
    tcm.getColumn(3).setPreferredWidth(170);

    JScrollPane scrollPane = new JScrollPane(batchQuantTable_);
    scrollPane.setPreferredSize(new Dimension(650, 200));
    quantifyingBatchPanel_.add(scrollPane,BorderLayout.CENTER);
    JPanel quantProgressPanel = new JPanel();
    quantProgressPanel.setLayout(new GridBagLayout());
    quantifyingBatchPanel_.add(quantProgressPanel,BorderLayout.SOUTH);
    
    quantifyingBatchLabel_ = new JLabel("Quantifying");
    quantifyingBatchLabel_.setToolTipText(TooltipTexts.QUANTITATION_STATUS_TEXT);
    quantProgressPanel.add(quantifyingBatchLabel_,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    progressBatchBar_ = new JProgressBar();
    progressBatchBar_.setMaximum(100);
    progressBatchBar_.setToolTipText(TooltipTexts.QUANTITATION_PROGRESS);
    quantProgressPanel.add(progressBatchBar_,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    Icon icon = new ImageIcon(LipidDataAnalyzer.class.getResource("/images/spinner.gif"));
    spinnerBatchLabel_ = new JLabel(icon);
    quantProgressPanel.add(spinnerBatchLabel_,new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    
    quantifyingBatchPanel_.setVisible(false);
  }
  
  private void createResultsMenu(){
    resultsMenu_ = new JPanel();
    resultsMenu_.setLayout(new GridBagLayout());
    jButtonResultFilesOpen_ = new JButton("Add Files", addFilesIcon_);
    jButtonResultFilesOpen_.addActionListener(this);
    jButtonResultFilesOpen_.setActionCommand("showResultsFilesChooser");
    jButtonResultFilesOpen_.setToolTipText(TooltipTexts.STATISTICS_ADD_RESULT_FILES);
    resultsMenu_.add(jButtonResultFilesOpen_,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    
    jButtonResultsDirOpen_ = new JButton("Add Dir",addFolderIcon_);
    jButtonResultsDirOpen_.addActionListener(this);
    jButtonResultsDirOpen_.setActionCommand("showResultsDirChooser");
    jButtonResultsDirOpen_.setToolTipText(TooltipTexts.STATISTICS_ADD_DIR);
    resultsMenu_.add(jButtonResultsDirOpen_,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    jButtonResultFilesRemove_ = new JButton("Remove",removeFilesIcon_);
    jButtonResultFilesRemove_.addActionListener(this);
    jButtonResultFilesRemove_.setActionCommand("removeResultFiles");
    jButtonResultFilesRemove_.setToolTipText(TooltipTexts.STATISTICS_REMOVE_SELECTION);
    resultsMenu_.add(jButtonResultFilesRemove_,new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));    
    
    jButtonResultFilesClean_ = new JButton("Remove all",removeFilesIcon_);
    jButtonResultFilesClean_.addActionListener(this);
    jButtonResultFilesClean_.setActionCommand("removeAllResultFiles");
    jButtonResultFilesClean_.setToolTipText(TooltipTexts.STATISTICS_REMOVE_ALL_SELECTION);
    resultsMenu_.add(jButtonResultFilesClean_,new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    jButtonResultAddToGroup_ = new JButton("Add to group",addToGroupIcon_);
    jButtonResultAddToGroup_.addActionListener(this);
    jButtonResultAddToGroup_.setActionCommand("addToGroup");
    jButtonResultAddToGroup_.setToolTipText(TooltipTexts.STATISTICS_ADD_TO_GROUP);
    resultsMenu_.add(jButtonResultAddToGroup_,new GridBagConstraints(4, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));

    analysisSelectionTablePanel_ = new JPanel();
    this.generateResultsAnalysisTablePane(new String[0][0]);
//    resultsMenu_.add(analysisSelectionTablePanel_);
    resultsMenu_.add(analysisSelectionTablePanel_,new GridBagConstraints(0, 1, 5, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    
    JPanel groupRtPanel = new JPanel();
    groupRtPanel.setLayout(new GridBagLayout());  
    JLabel sepRtLabel = new JLabel("Show hits with different RT separately: ");
    sepRtLabel.setToolTipText(TooltipTexts.STATISTICS_SEPARATE_RT);
    groupRtPanel.add(sepRtLabel,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    separateHitsByRT_  = new JCheckBox();
    separateHitsByRT_.setSelected(true);
    separateHitsByRT_.setActionCommand(CHANGE_SEPARATE_RT_STATUS);
    separateHitsByRT_.addActionListener(this);
    separateHitsByRT_.setToolTipText(TooltipTexts.STATISTICS_SEPARATE_RT);
    groupRtPanel.add(separateHitsByRT_,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    rtGroupingTime_ = new JTextField(3);
    rtGroupingTime_.setHorizontalAlignment(JTextField.RIGHT);
    rtGroupingTime_.setInputVerifier(new DoubleVerifier(true));
    rtGroupingTime_.setText("0.1");
    rtGroupingTime_.setToolTipText(TooltipTexts.STATISTICS_SEPARATE_RT);
    groupRtPanel.add(rtGroupingTime_,new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    rtTimeUnit_ = new JLabel("min");
    rtTimeUnit_.setToolTipText(TooltipTexts.STATISTICS_SEPARATE_RT);
    groupRtPanel.add(rtTimeUnit_,new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    resultsMenu_.add(groupRtPanel,new GridBagConstraints(0, 2, 5, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 15, 0), 0, 0));    
    
    JPanel correctOrderPanel = new JPanel();
    correctOrderPanel.setLayout(new GridBagLayout());  
    JLabel quantLabel = new JLabel("Quant file (for correct analyte order):");
    quantLabel.setToolTipText(TooltipTexts.STATISTICS_CORRECT_ORDER_FILE);
    correctOrderPanel.add(quantLabel,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    correctOrderFile_ = new JTextField(35);
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
    resultsMenu_.add(correctOrderPanel,new GridBagConstraints(0, 3, 5, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 15, 0), 0, 0));    
    
    groupsPanel_= new GroupsPanel();
    resultsMenu_.add(groupsPanel_ ,new GridBagConstraints(0, 4, 5, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    JLabel label = new JLabel("Internal-standard prefix: ");
    label.setToolTipText(TooltipTexts.STATISTICS_IS_PREFIX);
    resultsMenu_.add(label ,new GridBagConstraints(0, 5, 2, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    internalStandardSelection_ = new JTextField(3);
    internalStandardSelection_.setText(Settings.getInternalStandardDefaultInput());
    internalStandardSelection_.setToolTipText(TooltipTexts.STATISTICS_IS_PREFIX);
    resultsMenu_.add(internalStandardSelection_ ,new GridBagConstraints(1, 5, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    label = new JLabel("External-standard prefix: ");
    resultsMenu_.add(label ,new GridBagConstraints(0, 6, 2, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    label.setToolTipText(TooltipTexts.STATISTICS_ES_PREFIX);
    externalStandardSelection_ = new JTextField(3);
    externalStandardSelection_.setText(Settings.getExternalStandardDefaultInput());
    externalStandardSelection_.setToolTipText(TooltipTexts.STATISTICS_ES_PREFIX);
    resultsMenu_.add(externalStandardSelection_ ,new GridBagConstraints(1, 6, 1, 1, 0.0, 0.0
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
    resultsMenu_.add(cutoffPanel,new GridBagConstraints(3, 5, 2, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));

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
    resultsMenu_.add(absQuantSetPanel,new GridBagConstraints(3, 6, 2, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    cutoffSettingsPanel_ = new CutoffSettingsPanel(TooltipTexts.STATISTICS_ADD_CUTOFF_SETTINGS);
    resultsMenu_.add(cutoffSettingsPanel_ ,new GridBagConstraints(0, 7, 5, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));  
    quantSettingsPanel_ = new AbsoluteQuantSettingsPanel (groupsPanel_);
    resultsMenu_.add(quantSettingsPanel_ ,new GridBagConstraints(0, 8, 5, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));  
    jButtonResultFilesAccept_ = new JButton("Accept");
    jButtonResultFilesAccept_.addActionListener(this);
    jButtonResultFilesAccept_.setActionCommand("acceptSelectedResultFiles");
    jButtonResultFilesAccept_.setToolTipText(TooltipTexts.STATISTICS_ACCEPT_SELECTION);
    resultsMenu_.add(jButtonResultFilesAccept_,new GridBagConstraints(0, 9, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    
  }
  
  private void createSingleQuantMenu(){
    this.singleQuantMenu_ = new JPanel();
    this.singleQuantMenu_.setLayout(new GridBagLayout());
    JPanel selectionPanel = new JPanel();
    selectionPanel.setLayout(new GridBagLayout());
    singleQuantMenu_.add(selectionPanel,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));

    mzXMLLabel = new JLabel("Raw file: ");
    mzXMLLabel.setToolTipText(TooltipTexts.QUANTITATION_SINGLE_RAW_FILE);
    selectionPanel.add(mzXMLLabel,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    this.selectedMzxmlFile = new JTextField(62);
    selectedMzxmlFile.setToolTipText(TooltipTexts.QUANTITATION_SINGLE_RAW_FILE);
//    selectedFile.setEnabled(false);
    
    selectionPanel.add(selectedMzxmlFile,new GridBagConstraints(1, 0, 6, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    
    jButtonMxXMLOpen = new JButton("Select");
    jButtonMxXMLOpen.addActionListener(this);
    jButtonMxXMLOpen.setActionCommand("showMzxmlFileChooser");
    jButtonMxXMLOpen.setToolTipText(TooltipTexts.QUANTITATION_SINGLE_RAW_FILE);
    selectionPanel.add(jButtonMxXMLOpen,new GridBagConstraints(7, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    chromLabel = new JLabel("Quantitation: ");
    chromLabel.setToolTipText(TooltipTexts.QUANTITATION_SINGLE_MASS_LIST);
    this.selectedQuantFile = new JTextField(62);
    selectedQuantFile.setToolTipText(TooltipTexts.QUANTITATION_SINGLE_MASS_LIST);
    selectionPanel.add(selectedQuantFile,new GridBagConstraints(1, 1, 6, 1, 0.0, 0.0
      ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    selectionPanel.add(chromLabel,new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    jButtonQuantOpen = new JButton("Select");
    jButtonQuantOpen.addActionListener(this);
    jButtonQuantOpen.setActionCommand("showQuantFileChooser");
    jButtonQuantOpen.setToolTipText(TooltipTexts.QUANTITATION_SINGLE_MASS_LIST);
    selectionPanel.add(jButtonQuantOpen,new GridBagConstraints(7, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));

    
    JPanel settingsPanel = new JPanel ();  
    settingsPanel.setLayout(new GridBagLayout());
    singleQuantMenu_.add(settingsPanel,new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(6, 6, 0, 0), 0, 0));

//    JLabel mzTolLabel = new JLabel("m/z-Tolerance: ");
//    settingsPanel.add(mzTolLabel,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
//        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
//    singleMzTol_ = new JTextField(4);
//    singleMzTol_.setText("0.02");
//    singleMzTol_.setHorizontalAlignment(JTextField.RIGHT);
//    settingsPanel.add(singleMzTol_,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
//        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    JLabel timeMinusTolLabel = new JLabel("Time before tol.: ");
    timeMinusTolLabel.setToolTipText(TooltipTexts.QUANTITATION_RET_BEFORE);
    settingsPanel.add(timeMinusTolLabel,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    singleTimeMinusTol_ = new JTextField(4);
    singleTimeMinusTol_.setText("5");
    singleTimeMinusTol_.setHorizontalAlignment(JTextField.RIGHT);
    singleTimeMinusTol_.setToolTipText(TooltipTexts.QUANTITATION_RET_BEFORE);
    settingsPanel.add(singleTimeMinusTol_,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    JLabel timeMinusUnit = new JLabel("min");
    timeMinusUnit.setToolTipText(TooltipTexts.QUANTITATION_RET_BEFORE);
    settingsPanel.add(timeMinusUnit,new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    JLabel timePlusTolLabel = new JLabel("Time after tol.: ");
    timePlusTolLabel.setToolTipText(TooltipTexts.QUANTITATION_RET_AFTER);
    settingsPanel.add(timePlusTolLabel,new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    singleTimePlusTol_ = new JTextField(4);
    singleTimePlusTol_.setText("5");
    singleTimePlusTol_.setHorizontalAlignment(JTextField.RIGHT);
    singleTimePlusTol_.setToolTipText(TooltipTexts.QUANTITATION_RET_AFTER);
    settingsPanel.add(singleTimePlusTol_,new GridBagConstraints(4, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    JLabel timePlusUnit = new JLabel("min");
    timePlusUnit.setToolTipText(TooltipTexts.QUANTITATION_RET_AFTER);
    settingsPanel.add(timePlusUnit,new GridBagConstraints(5, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    JLabel cutoffLabel = new JLabel("Rel. base-peak cutoff: ");
    cutoffLabel.setToolTipText(TooltipTexts.QUANTITATION_CUTOFF);
    settingsPanel.add(cutoffLabel,new GridBagConstraints(6, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    singleCutoff_ = new JTextField(4);
    singleCutoff_.setText(LipidomicsConstants.getBasePeakDefaultCutoff());
    singleCutoff_.setHorizontalAlignment(JTextField.RIGHT);
    singleCutoff_.setToolTipText(TooltipTexts.QUANTITATION_CUTOFF);
    settingsPanel.add(singleCutoff_,new GridBagConstraints(7, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    JLabel cutoffUnit = new JLabel("\u2030");
    cutoffUnit.setToolTipText(TooltipTexts.QUANTITATION_CUTOFF);
    settingsPanel.add(cutoffUnit,new GridBagConstraints(8, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    JLabel rtShiftLabel = new JLabel("RT-shift: ");
    rtShiftLabel.setToolTipText(TooltipTexts.QUANTIFICATION_RET_SHIFT);
    settingsPanel.add(rtShiftLabel,new GridBagConstraints(9, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    singleRTShift_ = new JTextField(4);
    singleRTShift_.setText("0.0");
    singleRTShift_.setHorizontalAlignment(JTextField.RIGHT);
    singleRTShift_.setToolTipText(TooltipTexts.QUANTIFICATION_RET_SHIFT);
    settingsPanel.add(singleRTShift_,new GridBagConstraints(10, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    JLabel rtShiftUnit = new JLabel("min");
    rtShiftUnit.setToolTipText(TooltipTexts.QUANTIFICATION_RET_SHIFT);
    settingsPanel.add(rtShiftUnit,new GridBagConstraints(11, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));

    
    JPanel settingsPanel2 = new JPanel ();
    settingsPanel2.setLayout(new GridBagLayout());
    singleQuantMenu_.add(settingsPanel2,new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(6, 6, 0, 0), 0, 0));
    isoValidation_ = new JCheckBox();
    isoValidation_.setSelected(true);
    isoValidation_.addItemListener(new LipidomicsItemListener("endisableAmountOfIsotopes"));
    isoValidation_.setToolTipText(TooltipTexts.QUANTITATION_ISOTOPES);
    settingsPanel2.add(isoValidation_,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    JLabel allowIsotopicValidation = new JLabel("Isotopic quantitation of ");
    allowIsotopicValidation.setToolTipText(TooltipTexts.QUANTITATION_ISOTOPES);
    settingsPanel2.add(allowIsotopicValidation,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    amountOfIsotopes_ = new JTextField(4);
    amountOfIsotopes_.setText("2");
    amountOfIsotopes_.setHorizontalAlignment(JTextField.RIGHT);
    amountOfIsotopes_.setToolTipText(TooltipTexts.QUANTITATION_ISOTOPES);
    settingsPanel2.add(amountOfIsotopes_,new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    JLabel isotopePeaks = new JLabel("isotopes where");
    isotopePeaks.setToolTipText(TooltipTexts.QUANTITATION_ISOTOPES);
    settingsPanel2.add(isotopePeaks,new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    amountOfMatchingSearchIsotopes_ = new JTextField(4);
    amountOfMatchingSearchIsotopes_.setText("1");
    amountOfMatchingSearchIsotopes_.setHorizontalAlignment(JTextField.RIGHT);
    amountOfMatchingSearchIsotopes_.setToolTipText(TooltipTexts.QUANTITATION_ISOTOPES);
//    amountOfMatchingSearchIsotopes_.setEnabled(false);
    settingsPanel2.add(amountOfMatchingSearchIsotopes_,new GridBagConstraints(4, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    JLabel isotopeMatchPeaks = new JLabel("isotopic peak(s) have to match");
    isotopeMatchPeaks.setToolTipText(TooltipTexts.QUANTITATION_ISOTOPES);
    settingsPanel2.add(isotopeMatchPeaks,new GridBagConstraints(5, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));

    
    JPanel settingsPanel3 = new JPanel ();
    settingsPanel3.setLayout(new GridBagLayout());
    singleQuantMenu_.add(settingsPanel3,new GridBagConstraints(0, 3, 1, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(6, 6, 0, 0), 0, 0));
    searchUnknownTime_ = new JCheckBox();
    searchUnknownTime_.setSelected(true);
    searchUnknownTime_.setToolTipText(TooltipTexts.QUANTITATION_RET_UNKNOWN);
    settingsPanel3.add(searchUnknownTime_,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    JLabel findUnknownMolecules = new JLabel("Find molecules where retention time is unknown");
    findUnknownMolecules.setToolTipText(TooltipTexts.QUANTITATION_RET_UNKNOWN);
    settingsPanel3.add(findUnknownMolecules,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    JLabel maxProcessors = new JLabel("Processors for quantitation:");
    maxProcessors.setToolTipText(TooltipTexts.QUANTITATION_PROCESSORS);
    settingsPanel3.add(maxProcessors,new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    nrProcessors_ = new JTextField(3);
    nrProcessors_.setText(String.valueOf(this.getAmountOfProcessorsPreferred()));
    nrProcessors_.setHorizontalAlignment(JTextField.RIGHT);
    nrProcessors_.setToolTipText(TooltipTexts.QUANTITATION_PROCESSORS);
    nrProcessors_.setInputVerifier(new IntegerMaxVerifier(true,1,getMaxProcessors()));
    settingsPanel3.add(nrProcessors_,new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    
    int yPositionCounterSettingsPanels = 4;
    if (Settings.useAlex()){
      JPanel settingsPanel4 = new JPanel();
      settingsPanel4.setLayout(new GridBagLayout());
      singleQuantMenu_.add(settingsPanel4,new GridBagConstraints(0, yPositionCounterSettingsPanels, 1, 1, 0.0, 0.0
          ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(6, 6, 0, 0), 0, 0));
      JLabel ionModeLabel = new JLabel("Ion mode:");
      ionModeLabel.setToolTipText(TooltipTexts.QUANTITATION_ION_MODE);
      settingsPanel4.add(ionModeLabel,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
          ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
      ionMode_ = new JComboBox<String>();
      ionMode_.addItem("+");
      ionMode_.addItem("-");
      ionMode_.setToolTipText(TooltipTexts.QUANTITATION_ION_MODE);
      ionMode_.setFont(SELECT_FIELD_FONT);
      settingsPanel4.add(ionMode_,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
          ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
      JLabel ionModeInfo = new JLabel("(selection required for Alex123 target lists)");
      ionModeInfo.setToolTipText(TooltipTexts.QUANTITATION_ION_MODE);
      settingsPanel4.add(ionModeInfo,new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
      yPositionCounterSettingsPanels++;
    }

    startQuantification = new JButton("Start Quantitation");
    startQuantification.addActionListener(this);
    startQuantification.setActionCommand("startQuantification");
    startQuantification.setToolTipText(TooltipTexts.QUANTITATION_START);
    singleQuantMenu_.add(startQuantification,new GridBagConstraints(0, yPositionCounterSettingsPanels, 1, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(6, 6, 0, 0), 0, 0));
    yPositionCounterSettingsPanels++;
    
    quantifyingPanel_ = new JPanel ();  
    quantifyingPanel_.setLayout(new GridBagLayout());
    singleQuantMenu_.add(quantifyingPanel_,new GridBagConstraints(0, yPositionCounterSettingsPanels, 1, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(6, 6, 0, 0), 0, 0));
    yPositionCounterSettingsPanels++;
    quantifyingLabel_ = new JLabel("Quantifying");
    quantifyingLabel_.setToolTipText(TooltipTexts.QUANTITATION_STATUS_TEXT);
    quantifyingPanel_.add(quantifyingLabel_,new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    progressBar_ = new JProgressBar();
    progressBar_.setMaximum(100);
    progressBar_.setToolTipText(TooltipTexts.QUANTITATION_PROGRESS);
    quantifyingPanel_.add(progressBar_,new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    Icon icon = new ImageIcon(LipidDataAnalyzer.class.getResource("/images/spinner.gif"));
    spinnerLabel_ = new JLabel(icon);
    quantifyingPanel_.add(spinnerLabel_,new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
    

    quantifyingPanel_.setVisible(false);

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
         this.selectedMzxmlFile.setText(text);
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
         this.selectedMzxmlDirectory_.setText(text);
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
               if (suffix.equalsIgnoreCase("xls")||suffix.equalsIgnoreCase("xlsx")){
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
         if (jButtonResultAbsQuant_.getText().equalsIgnoreCase("Remove abs settings")||jButtonResultCutoff_.getText().equalsIgnoreCase("Remove cutoff settings"))
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
        if (resultFileCandidates.length>0 && jButtonResultAbsQuant_.getText().equalsIgnoreCase("Remove abs settings")||jButtonResultCutoff_.getText().equalsIgnoreCase("Remove cutoff settings"))
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
      int[] selectedColumns = resultFilesDisplayTable.getSelectedRows();
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
        if (jButtonResultAbsQuant_.getText().equalsIgnoreCase("Remove abs settings")||jButtonResultCutoff_.getText().equalsIgnoreCase("Remove cutoff settings"))
          subCommand = "removeSettings";
      }
    }
    if (this.resultFiles_!=null && command.equalsIgnoreCase("removeAllResultFiles")){
      this.resultFiles_ = new Vector<File>();
      this.groupsPanel_.removeAllGroups();
      this.updateAnalysisSelectionTable();
      if (jButtonResultAbsQuant_.getText().equalsIgnoreCase("Remove abs settings")||jButtonResultCutoff_.getText().equalsIgnoreCase("Remove cutoff settings"))
        subCommand = "removeSettings";
    }
    if (command.equalsIgnoreCase("addToGroup")){
      if (resultFiles_!=null){
        int[] selectedColumns = resultFilesDisplayTable.getSelectedRows();
        if (selectedColumns!=null && selectedColumns.length>0){
          Vector<File> selectedFiles = new Vector<File>();
          for (int selectedColumn : selectedColumns){
            selectedFiles.add(resultFiles_.get(selectedColumn));
          }
          resultFilesDisplayTable.getSelectionModel().clearSelection();
          int amountGroups = groupsPanel_.getGroups().size();
          InputDialog dlg = new InputDialog(new JFrame(), "Add to group", "Enter the group name", "Group "+String.valueOf(amountGroups+1));
          String groupName = dlg.getEnteredText();
          if (groupName!=null&&groupName.length()>0){
            groupName = groupName.trim();
            groupsPanel_.addGroup(groupName,selectedFiles,getWidth());
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
      this.acceptResultFiles();
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
         this.selectedQuantFile.setText(text);
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
         this.selectedQuantDir_.setText(text);
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
      if (selectedMzxmlDirectory_.getText()!=null&&selectedMzxmlDirectory_.getText().length()>0 &&
          selectedQuantDir_.getText()!=null&&selectedQuantDir_.getText().length()>0){
        canStartQuantification = true;
        if (this.isoBatchValidation_.isSelected()){
          if (this.amountOfBatchIsotopes_.getText()==null||this.amountOfBatchIsotopes_.getText().length()<1)
            canStartQuantification = false;
          else{
            try{
              Integer.parseInt(this.amountOfBatchIsotopes_.getText());
            }catch(NumberFormatException nfx){
              canStartQuantification = false;
            }
          }
          if (this.amountOfMatchingBatchSearchIsotopes_.getText()==null||this.amountOfMatchingBatchSearchIsotopes_.getText().length()<1)
            canStartQuantification = false;
          else{
            try{
              int isos = Integer.parseInt(this.amountOfMatchingBatchSearchIsotopes_.getText());
              if (isos<0)
                canStartQuantification = false;
            }catch(NumberFormatException nfx){
              canStartQuantification = false;
            }
          }
        }
        try{
          if (batchCutoff_.getText()!=null && batchCutoff_.getText().length()>0){
            cutoff = Float.parseFloat(batchCutoff_.getText().replaceAll(",", "."));
            LipidomicsConstants.getInstance().setRelativeMS1BasePeakCutoff(batchCutoff_.getText());
          }  
        }catch(NumberFormatException ex){new WarningMessage(new JFrame(), "Error", "The cutoff value must be float format!"); canStartQuantification=false;}
        try{
          if (batchRTShift_.getText()!=null && batchRTShift_.getText().length()>0)
            rtShift = Float.parseFloat(batchRTShift_.getText().replaceAll(",", "."));
        }catch(NumberFormatException ex){new WarningMessage(new JFrame(), "Error", "The RT-shift value must be float format!"); canStartQuantification=false;}

//        if (this.searchUnknownBatchTime_.isSelected()){
//        }
      }
      float minusTimeTol = 0f;
      float plusTimeTol = 0f;
      if (this.batchTimeMinusTol_.getText()!=null && this.batchTimeMinusTol_.getText().length()>0){
        try{ minusTimeTol = Float.parseFloat(this.batchTimeMinusTol_.getText());} catch (NumberFormatException nfx){new WarningMessage(new JFrame(), "Error", "The  \"Time before tol.\" value must be entered in float format!"); canStartQuantification=false;}
      }
      if (this.batchTimePlusTol_.getText()!=null && this.batchTimePlusTol_.getText().length()>0){
        try{ plusTimeTol = Float.parseFloat(this.batchTimePlusTol_.getText());} catch (NumberFormatException nfx){new WarningMessage(new JFrame(), "Error", "The  \"Time after tol.\" value must be entered in float format!"); canStartQuantification=false;}
      }
      if (canStartQuantification){      
        Vector<File> rawFiles = new Vector<File>();
        Vector<File> quantFiles = new Vector<File>();
        File rawDir = new File(this.selectedMzxmlDirectory_.getText());
        File quantDir = new File(this.selectedQuantDir_.getText());
        int amountOfIsotopes = 0;
        int isotopesMustMatch = 0;
        if (this.isoBatchValidation_.isSelected()){
          amountOfIsotopes = Integer.parseInt(this.amountOfBatchIsotopes_.getText());
          isotopesMustMatch = Integer.parseInt(this.amountOfMatchingBatchSearchIsotopes_.getText());
        }
//        if (this.searchUnknownBatchTime_.isSelected()){          
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
            if (this.ionModeBatch_!=null && ((String)ionModeBatch_.getSelectedItem()).equalsIgnoreCase("+"))
              ionMode = true;
            this.batchQuantTableModel_.clearFiles();
            this.batchQuantTableModel_.addFiles(pairs);
            this.quantifyingBatchLabel_.setText("Quantifying");
            this.progressBatchBar_.setValue(0);
            this.quantifyingBatchPanel_.setVisible(true);
            this.startBatchQuantification_.setEnabled(false);
            this.spinnerBatchLabel_.setVisible(true);
            batchQuantThread_ = new BatchQuantThread(this.batchQuantTable_, this.batchQuantTableModel_,this.progressBatchBar_, 
                this.quantifyingBatchLabel_,//Float.parseFloat(this.batchMzTol_.getText()),
                minusTimeTol,plusTimeTol,amountOfIsotopes,isotopesMustMatch,this.searchUnknownBatchTime_.isSelected(), cutoff, 
                rtShift, Integer.parseInt(nrProcessorsBatch_.getText()),ionMode,false);
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
        if (selectedMzxmlDirectory_.getText()==null||selectedMzxmlDirectory_.getText().length()<1){
          @SuppressWarnings("unused")
          WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "You must specify a raw, mzXML or chrom file directory");
//          wasWarningMessage=true;
        }else if (selectedQuantDir_.getText()==null||selectedQuantDir_.getText().length()<1){
          @SuppressWarnings("unused")
          WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "You must specify a quant file directory");
//          wasWarningMessage=true;
        } else if (this.isoBatchValidation_.isSelected()){
          if (this.amountOfBatchIsotopes_.getText()==null||this.amountOfBatchIsotopes_.getText().length()<1){
            @SuppressWarnings("unused")
            WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "If you select the isotope option, please specify the amount of isotopes (integer format)!");
//            wasWarningMessage=true;
          }else{
            try{
              Integer.parseInt(this.amountOfBatchIsotopes_.getText());
            }catch(NumberFormatException nfx){
              @SuppressWarnings("unused")
              WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "The number of isotopes must be integer format!");
//              wasWarningMessage=true;
            }
          }
        }
        if (this.amountOfMatchingBatchSearchIsotopes_.getText()==null||this.amountOfMatchingBatchSearchIsotopes_.getText().length()<0){
          @SuppressWarnings("unused")
          WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "If you select the find molecules without retention-time option, please specify the amount of matching isotopes (integer format)!");
//          wasWarningMessage=true;
        }else{
          try{
            int isotopes = Integer.parseInt(this.amountOfMatchingBatchSearchIsotopes_.getText());
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

//        if (!wasWarningMessage && this.searchUnknownBatchTime_.isSelected()){
//        }
      }
    }
    if (command.equalsIgnoreCase("startQuantification")){
      boolean canStartQuantification = false;
      float beforeTolerance = 0f;
      float afterTolerance = 0f;
      float cutoff = 0f;
      float rtShift = 0f;
      if (selectedMzxmlFile.getText()!=null)
        selectedMzxmlFile.setText(selectedMzxmlFile.getText().trim());
      if (selectedQuantFile.getText()!=null)
        selectedQuantFile.setText(selectedQuantFile.getText().trim());
      if (selectedMzxmlFile.getText()!=null&&selectedMzxmlFile.getText().length()>0 &&
          selectedQuantFile.getText()!=null&&selectedQuantFile.getText().length()>0){
        canStartQuantification = true;
        if (!StaticUtils.existsFile(selectedMzxmlFile.getText())){
          new WarningMessage(new JFrame(), "Error", "The raw file \""+StaticUtils.extractFileName(selectedMzxmlFile.getText())+"\" does not exist!"); 
          return;
        }else{
          if (selectedQuantFile.getText().length()>3){
            String suffix = selectedQuantFile.getText().substring(selectedQuantFile.getText().lastIndexOf("."));
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
          String suffix = selectedMzxmlFile.getText().substring(selectedMzxmlFile.getText().lastIndexOf(".")+1);
          if (!(suffix.equalsIgnoreCase("chrom")||suffix.equalsIgnoreCase("head")||suffix.equalsIgnoreCase("idx")||
              suffix.equalsIgnoreCase("rtt")||suffix.equalsIgnoreCase("raw")||suffix.equalsIgnoreCase("d")||suffix.equalsIgnoreCase("wiff")||
              suffix.equalsIgnoreCase(AbstractXMLSpectraReader.FILE_TYPE_MZ_XML)||
              suffix.equalsIgnoreCase(AbstractXMLSpectraReader.FILE_TYPE_MZ_ML))){
            new WarningMessage(new JFrame(), "Error", "For the raw data just files and directories with the suffix raw, mzXML, mzML, d, wiff, and chrom are allowed!"); 
            return;
          }
        }
        if (!StaticUtils.existsFile(selectedQuantFile.getText())){
          new WarningMessage(new JFrame(), "Error", "The mass list file \""+StaticUtils.extractFileName(selectedQuantFile.getText())+"\" does not exist!"); 
          canStartQuantification=false;
        }
        try{
          if (singleTimeMinusTol_.getText()!=null && singleTimeMinusTol_.getText().length()>0)
            beforeTolerance = Float.parseFloat(singleTimeMinusTol_.getText().replaceAll(",", "."));
        }catch(NumberFormatException ex){new WarningMessage(new JFrame(), "Error", "The  \"Time before tol.\" value must be entered in float format!"); canStartQuantification=false;}
        try{
          if (singleTimePlusTol_.getText()!=null && singleTimePlusTol_.getText().length()>0)
            afterTolerance = Float.parseFloat(singleTimePlusTol_.getText().replaceAll(",", "."));
        }catch(NumberFormatException ex){new WarningMessage(new JFrame(), "Error", "The  \"Time after tol.\" value must be entered in float format!"); canStartQuantification=false;}
        if (this.isoValidation_.isSelected()){
          if (this.amountOfIsotopes_.getText()==null||this.amountOfIsotopes_.getText().length()<1)
            canStartQuantification = false;
          else{
            try{
              Integer.parseInt(this.amountOfIsotopes_.getText());
            }catch(NumberFormatException nfx){
              canStartQuantification = false;
            }
          }
          if (this.amountOfMatchingSearchIsotopes_.getText()==null||this.amountOfMatchingSearchIsotopes_.getText().length()<1)
            canStartQuantification = false;
          else{
            try{
              int isos = Integer.parseInt(this.amountOfMatchingSearchIsotopes_.getText());
              if (isos<0)
                canStartQuantification = false;
            }catch(NumberFormatException nfx){
              canStartQuantification = false;
            }
          }
        }
        try{
          if (singleCutoff_.getText()!=null && singleCutoff_.getText().length()>0){
            cutoff = Float.parseFloat(singleCutoff_.getText().replaceAll(",", "."));
            LipidomicsConstants.getInstance().setRelativeMS1BasePeakCutoff(singleCutoff_.getText());
          }
        }catch(NumberFormatException ex){new WarningMessage(new JFrame(), "Error", "The cutoff value must be float format!"); canStartQuantification=false;}
        try{
          if (singleRTShift_.getText()!=null && singleRTShift_.getText().length()>0)
            rtShift = Float.parseFloat(singleRTShift_.getText().replaceAll(",", "."));
        }catch(NumberFormatException ex){new WarningMessage(new JFrame(), "Error", "The RT-shift value must be float format!"); canStartQuantification=false;}

//        if (this.searchUnknownTime_.isSelected()){
//        }
      }      
      if (canStartQuantification){
        this.quantifyingPanel_.setVisible(true);
        this.startQuantification.setEnabled(false);
        spinnerLabel_.setVisible(true);
        String fileToTranslate = selectedMzxmlFile.getText();
        readFromRaw_ = false;
        boolean threadStarted = false;
        boolean aborted = false;
        int amountOfIsotopes = 0;
        int isotopesMustMatch = 0;
        if (this.isoValidation_.isSelected()){
          amountOfIsotopes = Integer.parseInt(this.amountOfIsotopes_.getText());
          isotopesMustMatch = Integer.parseInt(this.amountOfMatchingSearchIsotopes_.getText());
        }
        
//        if (this.searchUnknownTime_.isSelected()){
//          isotopesMustMatch = Integer.parseInt(this.amountOfMatchingSearchIsotopes_.getText());
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
                this.progressBar_.setValue(5);
                this.quantifyingLabel_.setText("Translating to "+LipidomicsConstants.getIntermediateFileFormat());
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
            this.progressBar_.setValue(30);
            this.quantifyingLabel_.setText("Translating to chrom");
            mzToChromThread_ = new MzxmlToChromThread(fileToTranslate,Integer.parseInt(nrProcessors_.getText()));
            mzToChromThread_.start();
            threadStarted = true;
          }
        }
        if (!threadStarted && !aborted){
          boolean ionMode = false;
          if (this.ionMode_!=null && ((String)ionMode_.getSelectedItem()).equalsIgnoreCase("+"))
            ionMode = true;
          this.quantifyingLabel_.setText("Quantifying");
          this.progressBar_.setValue(70);
          quantThread_ = new QuantificationThread(selectedMzxmlFile.getText(), selectedQuantFile.getText(),LipidDataAnalyzer.getResultFilePath(selectedMzxmlFile.getText(), selectedQuantFile.getText()),
              beforeTolerance,afterTolerance,amountOfIsotopes,isotopesMustMatch,
              this.searchUnknownTime_.isSelected(),cutoff,rtShift,Integer.parseInt(nrProcessors_.getText()),ionMode,false);
          quantThread_.start();
          threadStarted = true;
        }else if (aborted){
          this.quantifyingPanel_.setVisible(false);
          this.startQuantification.setEnabled(true);
          spinnerLabel_.setVisible(false);          
        }
      }else{
//        boolean wasWarningMessage=false;
        if (selectedMzxmlFile.getText()==null||selectedMzxmlFile.getText().length()<1){
          @SuppressWarnings("unused")
          WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "You must specify a raw, mzXML or chrom file");
//          wasWarningMessage=true;
        }else if (selectedQuantFile.getText()==null||selectedQuantFile.getText().length()<1){
          @SuppressWarnings("unused")
          WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "You must specify a quant file");
//          wasWarningMessage=true;
        } else if (this.isoValidation_.isSelected()){
          if (this.amountOfIsotopes_.getText()==null||this.amountOfIsotopes_.getText().length()<1){
            @SuppressWarnings("unused")
            WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "If you select the isotope option, please specify the amount of isotopes that you want to quantify (integer format)!");
//            wasWarningMessage=true;
          }else{
            try{
              Integer.parseInt(this.amountOfIsotopes_.getText());
            }catch(NumberFormatException nfx){
              @SuppressWarnings("unused")
              WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "The number of isotopes must be integer format!");
//              wasWarningMessage=true;
            }
          }
          if (this.amountOfMatchingSearchIsotopes_.getText()==null||this.amountOfMatchingSearchIsotopes_.getText().length()<0){
            @SuppressWarnings("unused")
            WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "If you select the find molecules without retention-time option, please specify the amount of matching isotopes (integer format)!");
//            wasWarningMessage=true;
          }else{
            try{
              int isotopes = Integer.parseInt(this.amountOfMatchingSearchIsotopes_.getText());
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
//        if (!wasWarningMessage && this.searchUnknownTime_.isSelected()){
//        }
      }
    }
    if (command.equalsIgnoreCase("startDisplay")){
      if (selectedChromFile.getText()!=null)
        selectedChromFile.setText(selectedChromFile.getText().trim());
      if (selectedResultFile.getText()!=null)
        selectedResultFile.setText(selectedResultFile.getText().trim());
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
            this.readResultFile(selectedResultFile.getText());
            this.updateSheetSelectionList();
            this.params_ = null;
            RuleDefinitionInterface.clearCacheDir();
            lockRangeUpdateRequired_ = true;
          }
          catch (CgException e) {
            @SuppressWarnings("unused")
            WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "The chrom file you specified is not valid");
            e.printStackTrace();
          } catch (ExcelInputFileException e) {
            
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
    if (command.equalsIgnoreCase("updateQuantTolOfCurrentlySelected")){
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
      } else if (this.lockMzRange_.isSelected()){
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
            extractor.parseInput();
            jButtonResultAbsQuant_.setText("Remove abs settings");
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
            extractor.parseInput();
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
      rtGroupingTime_.setEnabled(separateHitsByRT_.isSelected());
      rtTimeUnit_.setEnabled(separateHitsByRT_.isSelected());
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
            spectrumPainter_.draw2DDiagram(svgGenerator);
            //this is for exporting chromatograms
            ////l2DPainter_.draw2DDiagram(svgGenerator);
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
      set.setRt(editRtDialog_.getRt());
      storeResultsToExcel();
      editRtDialog_ = null;
    }else if (command.equalsIgnoreCase("DeclineChangedRT")){
      editRtDialog_ = null;
    } else if (command.equalsIgnoreCase("AcceptExportSettings")){
      exportSettings_.setVisible(false);
      if (exportSettingsGroup_!=null)
        exportSettingsGroup_.setVisible(false);
    } 
    
    //TODO: remove when alternative isoLabel algo is finalized
//    else if (command.equalsIgnoreCase("AcceptOmegaExport")) {
//      Vector<IsotopicLabelVO> labelInfo = omegaExport_.getEnteredLabelInformation();
//      //TODO: this is only for debugging
//      for (IsotopicLabelVO label : labelInfo)
//        System.out.println(label);
//      
//      // get the directory where to store the omega mass list file
//      exportFileChooser_.setFileSelectionMode(JFileChooser.FILES_ONLY);
//      String fileName = "n-Masslist.xlsx";
//      String confirmDialogTitle = "\u03C9-RT export";
//      FileNameExtensionFilter filter = new FileNameExtensionFilter("Microsoft Office Excel Woorkbook (*.xlsx)","xlsx");
//      exportFileChooser_.setSelectedFile(new File(fileName));
//      exportFileChooser_.setFileFilter(filter);
//      if (JOptionPane.showConfirmDialog(this, "Only selected analytes with the current isotopes selected will be exported! Continue?",confirmDialogTitle,JOptionPane.YES_NO_OPTION) == JOptionPane.YES_OPTION){
//        int returnVal = exportFileChooser_.showSaveDialog(new JFrame());
//        if (returnVal != JFileChooser.APPROVE_OPTION)
//          return;
//        File fileToStore = exportFileChooser_.getSelectedFile();
//        @SuppressWarnings("rawtypes")
//        Vector results = HeatMapDrawing.checkFileStorage(fileToStore,"txt",resultsPanel_);
//        fileToStore = (File)results.get(0);
//        if ((Boolean)results.get(1))
//          exportOmegaMasslist(fileToStore,labelInfo);
//      }
//      exportFileChooser_.setSelectedFile(new File(""));
//    }
  }
  
  public void showExportSettingsDialog(boolean grouped){
    if (grouped)
      exportSettingsGroup_.setVisible(true);
    else
      exportSettings_.setVisible(true);
  }
  
  /**
   * updates and saves a LipidParameterSet object at a certain position in the displayed table
   * @param newOne the new LipidParameterSet object
   * @param position the position in the table where the object shall be updated
   */
  private void updateLipidParameterSet(LipidParameterSet newOne, int position){
    int originalPosition = resultPositionToOriginalLoopkup_.get(currentSelected_);
    result_.getIdentifications().get(currentSelectedSheet_).remove(originalPosition);
    result_.getIdentifications().get(currentSelectedSheet_).add(originalPosition,newOne);
    storeResultsToExcel();
  }
  
  /**
   * stores the current results back to the Excel file
   */
  private void storeResultsToExcel(){
    try {
      QuantificationResultExporter.writeResultsToExcel(selectedResultFile.getText(), result_);
      this.readResultFile(selectedResultFile.getText(),true);
      this.updateResultListSelectionTable();
      this.displayTable.changeSelection(this.currentSelected_, 1, false, false);
      listSelectionChanged(this.currentSelected_);
    } catch (Exception e) {
      e.printStackTrace();
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
       rangeColors = StaticUtils.createRangeColorVOs(param, ((LipidomicsTableModel)displayTable.getModel()).getMSnIdentificationName(currentSelected_),
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
  private void acceptResultFiles(){
    this.cleanupResultView();
    expDisplayNamesLookup_ = new Hashtable<String,String>();
    if (this.resultFiles_!=null&&this.resultFiles_.size()>0){
      AbsoluteSettingsVO absSettingVO = null;
      Hashtable<String,Double> cutoffValues = null;
      int maxCutoffIsotope = -1;
      if (jButtonResultAbsQuant_.getText().equalsIgnoreCase("Remove abs settings")){
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
      Hashtable<String,Vector<String>> correctAnalyteSequence = null;
      Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects = null;
      if (correctOrderFile_.getText()!=null && correctOrderFile_.getText().length()>0){
        try {
          boolean ionMode = false;
          if (this.ionModeOrder_!=null && ((String)ionModeOrder_.getSelectedItem()).equalsIgnoreCase("+"))
            ionMode = true;
          @SuppressWarnings("rawtypes")
          Vector quantInfo = QuantificationThread.getCorrectAnalyteSequence(correctOrderFile_.getText(),ionMode);
          classSequence = (LinkedHashMap<String,Integer>)quantInfo.get(0);
          correctAnalyteSequence = (Hashtable<String,Vector<String>>)quantInfo.get(1);
          quantObjects = (Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>>)quantInfo.get(3);
          
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
        analysisModule_.parseInput();
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
        e.printStackTrace();
        new WarningMessage(new JFrame(), "Error", "Some of the input files contain invalid information!");
      }
      catch (Exception e) {
        e.printStackTrace();
        new WarningMessage(new JFrame(), "Error", "The Excel returns the following failure: "+e.getMessage());
      }
        this.generateHeatMaps();
    }else{
      @SuppressWarnings("unused")
      WarningMessage dlg = new WarningMessage(new JFrame(), "Error", "Please specify the files to analyze! There have to be at least two!!");
    }
  }
  
  private void generateHeatMaps(){
    heatmaps_ = new Hashtable<String,HeatMapDrawing>();
    Hashtable<String,Hashtable<String,Hashtable<String,ResultCompVO>>> analysisResults = analysisModule_.getResults();
    Hashtable<String,Vector<String>> moleculeNames =  analysisModule_.getAllMoleculeNames();
    Hashtable<String,Hashtable<String,Integer>> corrTypeISLookup = analysisModule_.getCorrectionTypeISLookup();
    Hashtable<String,Hashtable<String,Integer>> corrTypeESLookup = analysisModule_.getCorrectionTypeESLookup();
    Vector<String> expNames = analysisModule_.getExpNamesInSequence();
    resultTabs_.removeAll();
    resultTabs_.addTab("Selection",resultsSelectionPanel_);
    resultTabs_.setToolTipTextAt(0, TooltipTexts.TABS_RESULTS_SELECTION);
    molBarCharts_ = new Hashtable<String,JTabbedPane>();
    
    // for the groupedValues
    Hashtable<String,Hashtable<String,Hashtable<String,ResultCompVO>>> groupedAnalysisResults = analysisModule_.getGroupedResults();
    groupHeatmaps_ = new Hashtable<String,HeatMapDrawing>();
    Hashtable<String,ResultDisplaySettings> displaySettingHash = new Hashtable<String,ResultDisplaySettings>();
    
    colorChooserDialog_ = new ColorChooserDialog(new JFrame(),"",expNames,groupsPanel_.getGroups(),this);
    int groupCount = 0;
    exportSettings_ = new ExportSettingsPanel(false,this);
    exportSettingsGroup_ = null;
    if (this.groupsPanel_.getGroups().size()>0)
      exportSettingsGroup_ = new ExportSettingsPanel(true,this);
    
    	//TODO: remove when alternative isoLabel algo is finalized
//    omegaExport_ = new OmegaExportDialog("\u03C9-RT export",analysisModule_.getIsoLabels(),this);
//    omegaExport_.setVisible(false);
    
    for (String molGroup : analysisResults.keySet()){
      JPanel aResultsViewPanel = new JPanel(new BorderLayout());
      JTabbedPane resultsViewTabs= new JTabbedPane();
      groupCount++;
      aResultsViewPanel.add(resultsViewTabs,BorderLayout.CENTER);
      Hashtable<String,Hashtable<String,ResultCompVO>> resultsOfOneGroup = analysisResults.get(molGroup);
      Vector<String> molNames = moleculeNames.get(molGroup);
      Hashtable<String,Integer> isLookup = new Hashtable<String,Integer> ();
      Hashtable<String,Integer> esLookup = new Hashtable<String,Integer> ();
      if (corrTypeISLookup.containsKey(molGroup))
        isLookup = corrTypeISLookup.get(molGroup);
      if (corrTypeESLookup.containsKey(molGroup))
        esLookup = corrTypeESLookup.get(molGroup);
      
      JPanel aPanel = new JPanel();
      aPanel.setLayout(new BorderLayout());
      boolean hasAbs = false;
      boolean hasProtein = false;
      boolean hasSampleWeight = false;
      boolean hasNeutralLipid = false;
      if (jButtonResultAbsQuant_.getText().equalsIgnoreCase("Remove abs settings")){
        hasAbs = true;
        try {
          if (this.quantSettingsPanel_.getSettingsVO().getVolumeSettings().size()>0 &&
              this.quantSettingsPanel_.getSettingsVO().getVolumeSettings().values().iterator().next().getProteinConc()!=null)
            hasProtein = true;
          if (this.quantSettingsPanel_.getSettingsVO().getVolumeSettings().size()>0 &&
              this.quantSettingsPanel_.getSettingsVO().getVolumeSettings().values().iterator().next().getNeutralLipidConc()!=null)
            hasNeutralLipid = true;
          if (this.quantSettingsPanel_.getSettingsVO().getVolumeSettings().size()>0 &&
              this.quantSettingsPanel_.getSettingsVO().getVolumeSettings().values().iterator().next().getSampleWeight()!=null)
            hasSampleWeight = true;

        }
        catch (AbsoluteSettingsInputException e) {
          // TODO Auto-generated catch block
          e.printStackTrace();
        }
      }
      ResultDisplaySettings displaySettings = new ResultDisplaySettings(analysisModule_.getISAvailability().get(molGroup),analysisModule_.getESAvailability().get(molGroup),isLookup,esLookup,hasAbs,
          hasSampleWeight,hasProtein,hasNeutralLipid);
      ResultSelectionSettings selectionSettings = new ResultSelectionSettings(null,molNames,true);
      ResultSelectionSettings combinedChartSettings = new ResultSelectionSettings(null,molNames,false);
      HeatMapDrawing drawing = new HeatMapDrawing(molGroup,resultsOfOneGroup, expNames,molNames, isLookup,esLookup,analysisModule_.getMaxIsotopesOfGroup(molGroup),analysisModule_.getModifications().get(molGroup), resultStatus_,this,molGroup,null,
          displaySettings,selectionSettings,combinedChartSettings,exportSettings_,analysisModule_.getRtTolerance());
      
      //TODO: remove when alternative isoLabel algo is finalized
//      HeatMapDrawing drawing = new HeatMapDrawing(molGroup,resultsOfOneGroup, expNames,molNames, isLookup,esLookup,analysisModule_.getMaxIsotopesOfGroup(molGroup),analysisModule_.getModifications().get(molGroup), resultStatus_,this,molGroup,null,
//          displaySettings,selectionSettings,combinedChartSettings,exportSettings_,omegaExport_, analysisModule_.getRtTolerance());
      displaySettings.addActionListener(drawing);
      selectionSettings.addActionListener(drawing);
      combinedChartSettings.addActionListener(drawing);
      JScrollPane scrollPane = new JScrollPane(drawing);
      int[] widthAndHeight = getScrollPaneWidthAndHeight(drawing);
      scrollPane.setPreferredSize(new Dimension(widthAndHeight[0], widthAndHeight[1]));
      heatmaps_.put(molGroup, drawing);
      if (expNames.size()>25)
        aPanel.add(scrollPane,BorderLayout.CENTER);
      else
        aPanel.add(scrollPane,BorderLayout.WEST);
      resultsViewTabs.addTab("Heatmap", aPanel);
      resultsViewTabs.setToolTipTextAt(0, TooltipTexts.TABS_RESULTS_HEATMAP+molGroup+"</html>");
      JPanel barChartPanel = new JPanel();
      resultsViewTabs.addTab("Bar-chart", barChartPanel);
      resultsViewTabs.setToolTipTextAt(1, TooltipTexts.TABS_RESULTS_BARCHART+molGroup+"</html>");
      if (this.groupsPanel_.getGroups().size()>0){
        Hashtable<String,Hashtable<String,ResultCompVO>> groupedResultsOfOneGroup = groupedAnalysisResults.get(molGroup);
        JPanel groupPanel = new JPanel();
        groupPanel.setLayout(new BorderLayout());
        HeatMapDrawing groupDrawing = new HeatMapDrawing(molGroup,groupedResultsOfOneGroup, this.groupsPanel_.getGroups(),molNames, isLookup,esLookup,analysisModule_.getMaxIsotopesOfGroup(molGroup),analysisModule_.getModifications().get(molGroup), resultStatus_,this,molGroup,
            drawing, displaySettings,selectionSettings,combinedChartSettings,exportSettingsGroup_,analysisModule_.getRtTolerance());
        //TODO: remove when alternative isoLabel algo is finalized
//        HeatMapDrawing groupDrawing = new HeatMapDrawing(molGroup,groupedResultsOfOneGroup, this.groupsPanel_.getGroups(),molNames, isLookup,esLookup,analysisModule_.getMaxIsotopesOfGroup(molGroup),analysisModule_.getModifications().get(molGroup), resultStatus_,this,molGroup,
//            drawing, displaySettings,selectionSettings,combinedChartSettings,exportSettingsGroup_,omegaExport_,analysisModule_.getRtTolerance());
        displaySettings.addActionListener(groupDrawing);
        selectionSettings.addActionListener(groupDrawing);
        combinedChartSettings.addActionListener(groupDrawing);
        JScrollPane groupScrollPane = new JScrollPane(groupDrawing);
        int[] groupWidthAndHeight = getScrollPaneWidthAndHeight(groupDrawing);
        scrollPane.setPreferredSize(new Dimension(groupWidthAndHeight[0], groupWidthAndHeight[1]));
        groupHeatmaps_.put(molGroup, groupDrawing);
        groupPanel.add(groupScrollPane,BorderLayout.WEST);
        resultsViewTabs.addTab("Group-Heatmap", groupPanel);
        resultsViewTabs.setToolTipTextAt(2, TooltipTexts.TABS_RESULTS_HEATMAP_GROUP+molGroup+"</html>");
        JPanel groupBarChartPanel = new JPanel();
        resultsViewTabs.addTab("Group bar-chart", groupBarChartPanel);
        resultsViewTabs.setToolTipTextAt(3, TooltipTexts.TABS_RESULTS_BARCHART_GROUP+molGroup+"</html>");
      }
      displaySettingHash.put(molGroup, displaySettings);
      molBarCharts_.put(molGroup, resultsViewTabs);
//        resultTabs_.addTab(molGroup,aPanel);
      resultTabs_.addTab(molGroup,aResultsViewPanel);
      resultTabs_.setToolTipTextAt(groupCount, TooltipTexts.TABS_RESULTS_GROUP+molGroup+"</html>");      
    }
    if (analysisResults.size()>0)
      resultTabs_.setSelectedIndex(1);
    if (analysisResults.size()>1){
      classOverviewPanel_ = new ClassesOverviewPanel(expNames,this,analysisModule_,displaySettingHash,heatmaps_,corrTypeISLookup,corrTypeESLookup,analysisResults,colorChooserDialog_);
      resultTabs_.addTab("Overview", classOverviewPanel_);
      resultTabs_.setToolTipTextAt(groupCount+1, TooltipTexts.TABS_RESULTS_OVERVIEW);
      if (this.groupsPanel_.getGroups().size()>0){
        classOverviewGroupPanel_ = new ClassesOverviewPanel(groupsPanel_.getGroups(),this,analysisModule_,displaySettingHash,heatmaps_,corrTypeISLookup,corrTypeESLookup,groupedAnalysisResults,colorChooserDialog_);
        resultTabs_.addTab("O.view-Group", classOverviewGroupPanel_);
        resultTabs_.setToolTipTextAt(groupCount+2, TooltipTexts.TABS_RESULTS_OVERVIEW_GROUPS);
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
      
      if (this.lockMzRange_.isSelected()) {
        startFloat = Float.parseFloat(displayMzStart_.getText());
        stopFloat = Float.parseFloat(displayMzStop_.getText());
        if (startFloat>stopFloat) {
          new WarningMessage(new JFrame(), "Error", "The stop value of the \""+DISPLAY_LOCK_MZ_TEXT+"\" cannot be smaller than the start m/z value");
          return;
        }
      } else {
        startFloat = currentIsotopicMass-Float.parseFloat(this.displayMinusTolerance_.getText());
        stopFloat = currentIsotopicMass+Float.parseFloat(this.displayPlusTolerance_.getText());
      }
      float startRt = 0f;
      if (this.displayRtStart_.getText()!=null && this.displayRtStart_.getText().length()>0)
        startRt = Float.parseFloat(this.displayRtStart_.getText());
      startRt = 60f*startRt;
      float stopRt = 0f;
      if (this.displayRtStop_.getText()!=null && this.displayRtStop_.getText().length()>0)
        stopRt = Float.parseFloat(this.displayRtStop_.getText());
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
      
      if (this.show2D_.isSelected()){
        java.awt.Panel l2DPainterPanel = new java.awt.Panel();
        Lipidomics2DPainter l2DPainter = null;
        if (params!=null){
          String[] rawLines2D = rawLines;
          float startFloat2D = startFloat;
          float stopFloat2D = stopFloat;
          if (this.lockMzRange_.isSelected() && ((currentIsotopicMass-2*params.LowerMzBand)<startFloat || stopFloat<(currentIsotopicMass+2*params.UpperMzBand))){
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

    frame_ = new MainFrame(new LipidDataAnalyzer(), 1050, 950);
    frame_.setTitle(getFrameTitleString());
  }
  
  private static String getFrameTitleString(){
    String titleString = "Lipid Data Analyzer "+Settings.VERSION+"   "+LipidomicsConstants.getCurrentMSMachine()+" settings ";
    String fragSelected = Settings.getFragmentSettingsString();
    if (fragSelected!=null) titleString += " "+fragSelected;    
    titleString += "         \u00A9 2023 - J\u00fcrgen Hartler, Andreas Ziegl, Gerhard G Thallinger, Leonida M Lamp - GNU GPL v3 license";
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
        this.startQuantification.setEnabled(true);
        spinnerLabel_.setVisible(false);
      }
      this.rawmzThread_ = null;
      this.progressBar_.setValue(30);
      this.quantifyingLabel_.setText("Translating to chrom");
      System.out.println("Translating to chrom from timer");
      this.readFromRaw_ = true;
      String filePath = selectedMzxmlFile.getText();
      String suffix = filePath.substring(filePath.lastIndexOf("."));
      if (suffix.equalsIgnoreCase(".wiff")){
        Vector<File> filesToTranslate = BatchQuantThread.getMzXMLFilesOfWiffConversion(filePath);
        if (filesToTranslate.size()==1){
          selectedMzxmlFile.setText(filesToTranslate.get(0).getAbsolutePath());
          mzToChromThread_ = new MzxmlToChromThread(filesToTranslate.get(0).getAbsolutePath(),Integer.parseInt(nrProcessors_.getText()));
          mzToChromThread_.start();          
        } else {
          Vector<RawQuantificationPairVO> pairs = new Vector<RawQuantificationPairVO>();
          File quantFile = new File(selectedQuantFile.getText());
          for (File rawFile: filesToTranslate){
            pairs.add(new RawQuantificationPairVO(rawFile,quantFile,true));
          }
          int amountOfIsotopes = 0;
          int isotopesMustMatch = 0;
          if (this.isoValidation_.isSelected()){
            amountOfIsotopes = Integer.parseInt(this.amountOfIsotopes_.getText());
            isotopesMustMatch = Integer.parseInt(this.amountOfMatchingSearchIsotopes_.getText());
          }
          boolean ok = true;
          float cutoff = 0f;
          float rtShift = 0f;
          try{
            if (singleCutoff_.getText()!=null && singleCutoff_.getText().length()>0)
              cutoff = Float.parseFloat(singleCutoff_.getText().replaceAll(",", "."));
          }catch(NumberFormatException ex){new WarningMessage(new JFrame(), "Error", "The cutoff value must be float format!"); ok=false;}
          try{
            if (singleRTShift_.getText()!=null && singleRTShift_.getText().length()>0)
              rtShift = Float.parseFloat(singleRTShift_.getText().replaceAll(",", "."));
          }catch(NumberFormatException ex){new WarningMessage(new JFrame(), "Error", "The RT-shift value must be float format!"); ok=false;}
          if (ok){
            boolean ionMode=false;
            if (this.ionMode_!=null && ((String)ionMode_.getSelectedItem()).equalsIgnoreCase("+"))
              ionMode = true;
            batchQuantTableModel_.clearFiles();
            batchQuantTableModel_.addFiles(pairs);
            batchQuantThread_ = new BatchQuantThread(this.batchQuantTable_, this.batchQuantTableModel_,this.progressBatchBar_, 
                this.quantifyingBatchLabel_,//Float.parseFloat(this.batchMzTol_.getText()),
                Float.parseFloat(this.singleTimeMinusTol_.getText()),Float.parseFloat(this.singleTimePlusTol_.getText()),
                amountOfIsotopes,isotopesMustMatch,this.searchUnknownTime_.isSelected(), cutoff, 
                rtShift, Integer.parseInt(nrProcessors_.getText()),ionMode,false);
            this.quantifyingBatchLabel_.setText("Quantifying");
            this.progressBatchBar_.setValue(0);
            this.quantifyingBatchPanel_.setVisible(true);
            this.startBatchQuantification_.setEnabled(false);
            this.spinnerBatchLabel_.setVisible(true);

            mainTabs.setSelectedIndex(1);
            batchQuantThread_.start();
          }
        }
      }else{
        String mzXMLFilePath = selectedMzxmlFile.getText().substring(0,selectedMzxmlFile.getText().lastIndexOf("."))+"."+LipidomicsConstants.getIntermediateFileFormat();
        mzToChromThread_ = new MzxmlToChromThread(mzXMLFilePath,Integer.parseInt(nrProcessors_.getText()));
        mzToChromThread_.start();
      }
    }
    if (this.mzToChromThread_!=null && this.mzToChromThread_.finished()){
      if (this.mzToChromThread_.getErrorString()!=null&&this.mzToChromThread_.getErrorString().length()>0){
        @SuppressWarnings("unused")
        WarningMessage dlg = new WarningMessage(new JFrame(), "Error", mzToChromThread_.getErrorString());
        this.startQuantification.setEnabled(true);
        spinnerLabel_.setVisible(false);
      }
      if (readFromRaw_){
        this.readFromRaw_ = false;
        RawToMzxmlThread.deleteMzXMLFiles(selectedMzxmlFile.getText().substring(0,selectedMzxmlFile.getText().lastIndexOf("."))+"."+LipidomicsConstants.getIntermediateFileFormat());
      }

      boolean ok = true;
      if (mzToChromThread_.isPolaritySwitched()){
        if (selectedQuantFile.getText().indexOf(GlobalConstants.CHROMATOGRAM_HEADER_FILE_POLARITY_POSITIVE)!=-1){
          String newRawName = selectedMzxmlFile.getText();
          newRawName = newRawName.substring(0,newRawName.lastIndexOf("."))+RawToChromThread.FILE_SUFFIX_POLARITY_POSITIVE+".chrom";
          selectedMzxmlFile.setText(newRawName);
        } else if (selectedQuantFile.getText().indexOf(GlobalConstants.CHROMATOGRAM_HEADER_FILE_POLARITY_NEGATIVE)!=-1){
          String newRawName = selectedMzxmlFile.getText();
          newRawName = newRawName.substring(0,newRawName.lastIndexOf("."))+RawToChromThread.FILE_SUFFIX_POLARITY_NEGATIVE+".chrom";
          selectedMzxmlFile.setText(newRawName);
        } else {
          ok = false;
          new WarningMessage(new JFrame(), "Error", "This is polarity switched data! It is not clear which of the two generated chrom files shall be quantified! Please use "+GlobalConstants.CHROMATOGRAM_HEADER_FILE_POLARITY_POSITIVE+" or "+GlobalConstants.CHROMATOGRAM_HEADER_FILE_POLARITY_NEGATIVE+" in the file name!");
          this.quantifyingPanel_.setVisible(false);
        }
      }
      this.mzToChromThread_ = null;
      this.quantifyingLabel_.setText("Quantifying");
      System.out.println("Quantifying from thread");
      this.progressBar_.setValue(70);
      int amountOfIsotopes = 0;
      int isotopesMustMatch = 0;
      if (this.isoValidation_.isSelected()){
        amountOfIsotopes = Integer.parseInt(this.amountOfIsotopes_.getText());
        isotopesMustMatch = Integer.parseInt(this.amountOfMatchingSearchIsotopes_.getText());
      }      
      float cutoff = 0f;
      float rtShift = 0f;
      try{
        if (singleCutoff_.getText()!=null && singleCutoff_.getText().length()>0)
          cutoff = Float.parseFloat(singleCutoff_.getText().replaceAll(",", "."));
      }catch(NumberFormatException ex){new WarningMessage(new JFrame(), "Error", "The cutoff value must be float format!"); ok=false;}
      try{
        if (singleRTShift_.getText()!=null && singleRTShift_.getText().length()>0)
          rtShift = Float.parseFloat(singleRTShift_.getText().replaceAll(",", "."));
      }catch(NumberFormatException ex){new WarningMessage(new JFrame(), "Error", "The RT-shift value must be float format!"); ok=false;}
      if (ok){
        boolean ionMode = false;
        if (this.ionMode_!=null && ((String)ionMode_.getSelectedItem()).equalsIgnoreCase("+"))
          ionMode = true;
        quantThread_ = new QuantificationThread(selectedMzxmlFile.getText(), selectedQuantFile.getText(),LipidDataAnalyzer.getResultFilePath(selectedMzxmlFile.getText(), selectedQuantFile.getText()),// Float.parseFloat(this.singleMzTol_.getText()),
            Float.parseFloat(this.singleTimeMinusTol_.getText()),Float.parseFloat(this.singleTimePlusTol_.getText()),amountOfIsotopes,isotopesMustMatch,this.searchUnknownTime_.isSelected(),
            cutoff,rtShift,Integer.parseInt(this.nrProcessors_.getText()),ionMode,false);
        quantThread_.start();
      }else{
        this.startQuantification.setEnabled(true);
        spinnerLabel_.setVisible(false);
      }
    }
    if (this.quantThread_!=null && !this.quantThread_.finished() && 
        (this.quantThread_.getErrorString()==null||this.quantThread_.getErrorString().length()==0)){
      if (quantThread_.getTotalAmountOfLipids()>0&&quantThread_.getCurrentLipidCount()>0){
        this.progressBar_.setValue(70+((30*(quantThread_.getCurrentLipidCount()-1))/quantThread_.getTotalAmountOfLipids()));
        quantifyingLabel_.setText("Quantifying "+quantThread_.getCurrentLipid()+" ("+quantThread_.getCurrentLipidCount()+"/"+quantThread_.getTotalAmountOfLipids()+")");
      }  
    }
    if (this.quantThread_!=null && this.quantThread_.finished()){
      if (this.quantThread_.getErrorString()!=null&&this.quantThread_.getErrorString().length()>0){
        @SuppressWarnings("unused")
        WarningMessage dlg = new WarningMessage(new JFrame(), "Error", quantThread_.getErrorString());
        this.startQuantification.setEnabled(true);
        spinnerLabel_.setVisible(false);
      }
      this.quantThread_ = null;
      this.progressBar_.setValue(100);
      quantifyingLabel_.setText("Finished");
      this.startQuantification.setEnabled(true);
      spinnerLabel_.setVisible(false);
    }
    if (this.batchQuantThread_!=null && this.batchQuantThread_.finished()){
      this.batchQuantThread_ = null;
      this.progressBatchBar_.setValue(100);
      quantifyingBatchLabel_.setText("Finished");
      this.startBatchQuantification_.setEnabled(true);
      spinnerBatchLabel_.setVisible(false);
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
    LipidomicsTableModel model = new LipidomicsTableModel(lipidsOrdered,lipids,showMSnNames_.isSelected(),resultsShowModification_.get(selectedSheet_.getSelectedItem()));
    resultPositionToOriginalLoopkup_ = model.getPositionToOriginal(); 
    displayTable = new LipidomicsJTable(model, new LipidomicsTableCellRenderer(),
        reader_.getHighestMsLevel()>1&&reader_.getMsmsType().equalsIgnoreCase(ChromatogramReader.CHROMATOGRAM_HEADER_FILE_MSMS_TYPE_PRECURSOR), orderType,
        QuantificationThread.hasRtInfo(result_.getIdentifications()),this);
    listSelectionModel = displayTable.getSelectionModel();
    displayTable.setSelectionModel(listSelectionModel);
    listSelectionModel.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
    tablePane = new JScrollPane(displayTable);
    tablePane.setPreferredSize(new Dimension(420, 130));
    tablePanel_.add(tablePane,BorderLayout.CENTER);
    tablePanel_.invalidate();
    tablePanel_.updateUI();
    selectionPane.setVisible(true);    
  }
  
  private void updateAnalysisSelectionTable(){
    analysisSelectionTablePanel_.remove(analysisTablePane);
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
    int columnWidth = 950;
    int tableWidth = columnWidth-3;
    if (tableData.length>24)
      tableWidth = columnWidth-18;
    String[] columnNames = { "file name", "directory" };
    resultFilesDisplayTable = new JTable(tableData, columnNames);
    resultListSelectionModel = resultFilesDisplayTable.getSelectionModel();
    resultFilesDisplayTable.setSelectionModel(resultListSelectionModel);
    resultListSelectionModel.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);  
    analysisTablePane = new JScrollPane(resultFilesDisplayTable);
    analysisTablePane.setPreferredSize(new Dimension(columnWidth, 300));
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
      int columnWidthOne = (int)((double)tableWidth*percentOne);
      int columnWidthTwo = tableWidth-columnWidthOne;
      resultFilesDisplayTable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
      resultFilesDisplayTable.getColumnModel().getColumn(0).setPreferredWidth(columnWidthOne);
      resultFilesDisplayTable.getColumnModel().getColumn(1).setPreferredWidth(columnWidthTwo);
    }
    analysisSelectionTablePanel_.add(analysisTablePane);
    
  }
   
  
  private void change2DTypeState(ItemEvent e, String command){
    if (command.equalsIgnoreCase("DisplayModeRaw") && this.l2DPainter_!=null){
      this.l2DPainter_.setRaw(true);
      this.l2DPainter_.repaint();
    }
    if (command.equalsIgnoreCase("DisplayModeSmooth") && this.l2DPainter_!=null){
      this.l2DPainter_.setRaw(false);
      this.l2DPainter_.repaint();
    }
    if (command.equalsIgnoreCase("show2dChanged")){
      if (this.params_!=null){
        if (this.show2D_.isSelected())
          this.initMS1OrMS2View(params_);
        else
          this.initANewViewer(params_);
      }
///        this.initANewViewer(params_);
    }
    if (command.equalsIgnoreCase("DisplayModeAbundance")){
      this.spectrumPainter_.setRelativeValues(relAbund_.isSelected());
      this.spectrumPainter_.repaint();
    }
    if (command.equalsIgnoreCase("endisableAmountOfIsotopes")){
      if (this.isoValidation_.isSelected()){
        this.amountOfIsotopes_.setEnabled(true);
        this.amountOfMatchingSearchIsotopes_.setEnabled(true);
      }else{
        this.amountOfIsotopes_.setEnabled(false);
        this.amountOfMatchingSearchIsotopes_.setEnabled(false);
      }  
    }
    if (command.equalsIgnoreCase("endisableAmountOfBatchIsotopes")){
      if (this.isoBatchValidation_.isSelected()){
        this.amountOfBatchIsotopes_.setEnabled(true);
        this.amountOfMatchingBatchSearchIsotopes_.setEnabled(true);
      }else{
        this.amountOfBatchIsotopes_.setEnabled(false);
        this.amountOfMatchingBatchSearchIsotopes_.setEnabled(false);
      }  
    }
//    if (command.equalsIgnoreCase("searchForUnknownRetentionTime")){
//      if (this.searchUnknownTime_.isSelected())
//        this.amountOfMatchingSearchIsotopes_.setEnabled(true);
//      else
//        this.amountOfMatchingSearchIsotopes_.setEnabled(false);
//    }
//    if (command.equalsIgnoreCase("searchForUnknownBatchRetentionTime")){
//      if (this.searchUnknownBatchTime_.isSelected())
//        this.amountOfMatchingBatchSearchIsotopes_.setEnabled(true);
//      else
//        this.amountOfMatchingBatchSearchIsotopes_.setEnabled(false);
//    }   
    if (command.equalsIgnoreCase("ChangeIsotope") && !this.displaysMs2_){
      if (e.getStateChange()==ItemEvent.SELECTED){
        if (this.params_!=null && l2DPainter_!=null){
          if (!this.lockMzRange_.isSelected() || this.viewer_==null || lockRangeUpdateRequired_){
            this.initANewViewer(params_,this.l2DPainter_.getAllSelectedProbes());
            lockRangeUpdateRequired_ = false;
          }else {
            updateViewForLockMz(params_);
          }        
        }
      }
    } else if (command.equalsIgnoreCase("showMSnNames")){
      updateResultListSelectionTable();
    } else if (command.equalsIgnoreCase(CHANGE_LOCK_MZ_RANGE)) {
      changeLockMzRange();
    }
    
  }
  
  private class LipidomicsItemListener implements java.awt.event.ItemListener
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
        this.displayTable.changeSelection(this.currentSelected_, 1, false, false);
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
    for (String sheetName : result_.getIdentifications().keySet()){ 
      selectedSheet_.addItem(sheetName);
    }
    selectedSheet_.setToolTipText(TooltipTexts.DISPLAY_SELECT_CLASS);
    tablePanel_.add(selectedSheet_,BorderLayout.NORTH);
    tablePanel_.invalidate();
    tablePanel_.updateUI();
    selectionPane.setVisible(true);
  }

  public boolean heatMapClicked(String experimentName, String resultFilePath,  String moleculeNameIn)
  {
    File resultsFile = new File (resultFilePath);
    String moleculeName = new String(moleculeNameIn);
    if (resultsFile.exists()&&resultsFile.isFile()){
      selectedResultFile.setText(resultFilePath);
      String chromFileBase = StaticUtils.extractChromBaseName(resultFilePath,experimentName);
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
        startDisplay.doClick();
        String sheetToSelect = resultTabs_.getTitleAt(resultTabs_.getSelectedIndex());
        selectedSheet_.setSelectedItem(sheetToSelect);
        String moelculeTableName = null;
//        String rtInTableName = null;
        int selection = -1;
        String[] molRtAndMod = StaticUtils.extractMoleculeRtAndModFromMoleculeName(moleculeName);
        moleculeName = molRtAndMod[0];
        String rt = molRtAndMod[1];
        for (int i=0;i!=this.displayTable.getRowCount();i++){
          String moleculeInTable = (String)this.displayTable.getSumLipidNameAt(i);
          if (moleculeInTable.startsWith(moleculeName)){
            boolean found = false;
            if (rt==null) found = true;
            else{
              String rtInTableString = moleculeInTable.substring(moleculeName.length()+1);
              if (rtInTableString.indexOf("_")!=-1) rtInTableString = rtInTableString.substring(0,rtInTableString.indexOf("_"));
              try{
                if (analysisModule_.isWithinRtGroupingBoundaries(Double.valueOf(rtInTableString), Double.valueOf(rt))){
                  found=true;
//                  rtInTableName = rtInTableString;
                }
              }catch(NumberFormatException nfx){}
            }
            if (found){
              //if show MS2 spectra is selected, try to find an adequate matching hit where MS2 spectra are present
              if (this.displaysMs2_){
                int j=i;
                boolean foundMsn = false;
                while (j<displayTable.getRowCount() && ((String)this.displayTable.getSumLipidNameAt(j)).startsWith(moleculeName)){
                  if (((LipidomicsTableModel)displayTable.getModel()).hasMS2Evidence(j)){
                    String moleculeInTableMsn = (String)this.displayTable.getSumLipidNameAt(j);
                    if (rt==null) foundMsn = true;
                    else{
                      String rtInTableString = moleculeInTableMsn.substring(moleculeName.length()+1);
                      if (rtInTableString.indexOf("_")!=-1) rtInTableString = rtInTableString.substring(0,rtInTableString.indexOf("_"));
                      try{
                        if (analysisModule_.isWithinRtGroupingBoundaries(Double.valueOf(rtInTableString), Double.valueOf(rt))){
                          foundMsn=true;
                        }
                      }catch(NumberFormatException nfx){}
                    }
                    if (foundMsn)
                      break;
                  }
                  j++;
                }
                if (foundMsn && i!=j){
                  i = j;
                  moleculeInTable = (String)this.displayTable.getSumLipidNameAt(i);
                }
              }
              moelculeTableName = moleculeInTable;
              selection = i;
              break;
            }
          }
        }
        if (moelculeTableName!=null){
          mainTabs.setSelectedIndex(3);
//          try {
//            Thread.sleep(100);
//          }
//          catch (InterruptedException e) {
//            // TODO Auto-generated catch block
//            e.printStackTrace();
//          }
//          ListSelectionEvent event2 = new ListSelectionEvent(displayTable,selection,selection+1,false);
//          displayTable.addSelectionInterval(selection,selection);
//          System.out.println("Handlers: "+handler.length);
//          .valueChanged(event2);
          this.makeDisplayRemoveOperations();
          this.displayTable.changeSelection(selection, 1, false, false);
          listSelectionChanged(selection);
        }else
          return false;
      }else{
        new WarningMessage(new JFrame(),"ERROR","The chrom file \""+experimentName+".chrom\" is not there!");
        return false;
      }  
    }else{
      new WarningMessage(new JFrame(),"ERROR","The result file \""+resultFilePath+"\" does not exist!");
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
        analysisModule_.getCorrectionTypeISLookup().get(groupName), analysisModule_.getCorrectionTypeESLookup().get(groupName),analysisModule_.getModifications().get(groupName),colorChooserDialog_),null,1);
    molBarCharts_.get(groupName).setSelectedIndex(1);
    return true; 
  }
  
  public boolean analyteGroupClicked(String moleculeName, String groupName, int maxIsotope, boolean rtGrouped, ResultDisplaySettingsVO settingVO, String prefUnit, String unit){
    molBarCharts_.get(groupName).remove(3);
    Hashtable<String,ResultCompVO> analysisResults = analysisModule_.getGroupedResults().get(groupName).get(moleculeName);
    Vector<String> groupNames = new Vector<String>(analysisModule_.getGroupNames());
    molBarCharts_.get(groupName).insertTab("Group bar-chart", null, new BarChartPainter(BarChartPainter.TYPE_MOLECULE,groupName,moleculeName,analysisResults,groupNames,this,false,true,
        getMaxIsotopeForSetting(settingVO, maxIsotope,analysisResults),rtGrouped, true, settingVO, prefUnit, unit,analysisModule_.getCorrectionTypeISLookup().get(groupName),
        analysisModule_.getCorrectionTypeESLookup().get(groupName),analysisModule_.getModifications().get(groupName),colorChooserDialog_)
    ,null,3);
    molBarCharts_.get(groupName).setSelectedIndex(3);
    return true; 
  }
  
  public static int getMaxIsotopeForSetting(ResultDisplaySettingsVO settingVO, int maxIsotope,Hashtable<String,ResultCompVO> analysisResults){
    int maxAppliIsotope = maxIsotope+1;
    if (!settingVO.getType().equalsIgnoreCase("relative to measured class amount") && 
        !settingVO.getType().equalsIgnoreCase("relative to total amount"))
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
        analysisModule_.getCorrectionTypeESLookup().get(groupName),analysisModule_.getModifications().get(groupName),colorChooserDialog_)
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
        analysisModule_.getCorrectionTypeESLookup().get(groupName),analysisModule_.getModifications().get(groupName),colorChooserDialog_)
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
        analysisModule_.getCorrectionTypeESLookup().get(groupName),analysisModule_.getModifications().get(groupName),colorChooserDialog_)
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
        analysisModule_.getCorrectionTypeESLookup().get(groupName),analysisModule_.getModifications().get(groupName),colorChooserDialog_)
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
                if (i==0)remainingParam.setRt(Calculator.FormatNumberToString(newIsotopicProbes.get(0).get(0).Peak/60f,2));
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
        acceptResultFiles();
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
                  String rt = "";
                  for (int k=0;k!=probes.size();k++){
                    Vector<CgProbe> isoProbes = probes.get(k);
                    Vector<Double> rts = new Vector<Double>();
                    for (CgProbe probe : isoProbes){
                      totalArea+=probe.Area;
                      rts.add((double)probe.Peak);
                    }
                    if (k==0 && hasRtInfo) rt = Calculator.FormatNumberToString(Calculator.mean(rts)/60d,2d);
                  }
                  if (probes.size()==0 && hasRtInfo) rt = templateParam.getRt();
// INFO: Settings.emptyEntriesForQuantAnalNotFound() is responsible for creating empty entries if the analyte cannot be quantified                  
                  if (totalArea>0 || Settings.emptyEntriesForQuantAnalNotFound()){
                    LipidParameterSet paramToQuantify = new LipidParameterSet(templateParam.Mz[0], templateParam.getName(),
                      templateParam.getDoubleBonds(), templateParam.getOhNumber(), templateParam.getModificationName(), rt, templateParam.getAnalyteFormula(),
                      templateParam.getModificationFormula(),templateParam.getCharge());
                      paramToQuantify.LowerMzBand = templateParam.LowerMzBand;
                    paramToQuantify.UpperMzBand = templateParam.UpperMzBand;
                    paramToQuantify.Area = totalArea;
                    paramToQuantify.setIsotopicProbes(probes);
                    String[] paramNameAndRt = StaticUtils.extractMoleculeRtAndModFromMoleculeName(paramToQuantify.getNameString());
                    boolean existsSameOne = false;
                    if (rt.length()>0){
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
                    if (rt.length()>0)
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
      acceptResultFiles();
      for (int i=0; i!=resultTabs_.getTabCount();i++){
        if (resultTabs_.getTitleAt(i).equalsIgnoreCase(groupName))
          resultTabs_.setSelectedIndex(i);
      }
    } catch (ExcelInputFileException eif){
      //Comment: the graphical Warning message is shown in the readResultFile itself
    }
  }
  

  public void eliminateAnalyteEverywhere(String groupName, Hashtable<String,String> selectedAnalytes, Vector<String> selectedMods, Vector<String> foundUpdateables){
    for (String updateablePath : foundUpdateables){
      Hashtable<String,String> modHash = new Hashtable<String,String>();
      for (String modName : selectedMods)modHash.put(modName, modName);
      Hashtable<String,Boolean> showMods = new Hashtable<String,Boolean>();
      Hashtable<String,Hashtable<String,String>> nameRtHash = new Hashtable<String,Hashtable<String,String>>();
      for (String selected : selectedAnalytes.keySet()){
        String[] nameAndRt = StaticUtils.extractMoleculeRtAndModFromMoleculeName(selected);
        if (nameAndRt[1]==null) nameAndRt[1]="";
        Hashtable<String,String> rts = new Hashtable<String,String>();
        if (nameRtHash.containsKey(nameAndRt[0])) rts = nameRtHash.get(nameAndRt[0]);
        rts.put(nameAndRt[1], nameAndRt[1]);
        nameRtHash.put(nameAndRt[0], rts);
      }
      try{
        QuantificationResult result = LDAResultReader.readResultFile(updateablePath, showMods);
        Vector<LipidParameterSet> oldParams = result.getIdentifications().get(groupName);
        Vector<LipidParameterSet> newParams = new Vector<LipidParameterSet>();
        for (LipidParameterSet param : oldParams){
          String[] nameAndRt = StaticUtils.extractMoleculeRtAndModFromMoleculeName(param.getNameString());
          if (nameAndRt[1]==null) nameAndRt[1]="";
          boolean remove = false;
          if (nameRtHash.containsKey(nameAndRt[0])){
            Hashtable<String,String> rts = nameRtHash.get(nameAndRt[0]);
            for (String rt : rts.keySet()){
              if (analysisModule_.neglectRtInformation(nameAndRt[0]) || analysisModule_.isWithinRtGroupingBoundaries(Double.valueOf(nameAndRt[1]), Double.valueOf(rt))){
                if (modHash.containsKey(param.getModificationName()))remove = true;
              }
            }
          }
          if (!remove) newParams.add(param);
        }
        result.getIdentifications().put(groupName, newParams);
        try {
          QuantificationResultExporter.writeResultsToExcel(updateablePath, result);
        }
        catch (Exception e) {
          e.printStackTrace();
        }
      } catch (ExcelInputFileException eif){
        //Comment: the graphical Warning message is shown in the readResultFile itself
      }
    }
    acceptResultFiles();
    for (int i=0; i!=resultTabs_.getTabCount();i++){
      if (resultTabs_.getTitleAt(i).equalsIgnoreCase(groupName))
        resultTabs_.setSelectedIndex(i);
    }    
  }
  
  public void addAnalyte(int position, AddAnalyteVO analyteDescrVO){
    Float exactMass = new Float(analyteDescrVO.getExactMass());
    LipidParameterSet set = new LipidParameterSet(exactMass, analyteDescrVO.getName(), 
        analyteDescrVO.getDoubleBonds(), analyteDescrVO.getOh(), analyteDescrVO.getModName(),
        analyteDescrVO.getRt(), analyteDescrVO.getFormula(), analyteDescrVO.getModFormula(),
        new Integer(analyteDescrVO.getCharge()) );
      set.LowerMzBand = LipidomicsConstants.getCoarseChromMzTolerance(exactMass);
      set.UpperMzBand = LipidomicsConstants.getCoarseChromMzTolerance(exactMass);
      set.Area = 0;
      result_.getIdentifications().get(this.currentSelectedSheet_).add(position, set);
      try {
       QuantificationResultExporter.writeResultsToExcel(selectedResultFile.getText(), result_);
       readResultFile(selectedResultFile.getText(),true);
       updateResultListSelectionTable();
       // this is that the selection gets updated in any case
       currentSelected_ = -1;
       displayTable.changeSelection(position, 1, false, false);
       listSelectionChanged(position);
    } catch (Exception e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
  }
  
  public void removeAnalyte(int[] indices){
    List<Integer> inds = new ArrayList<Integer>();
    for (int ind:indices)inds.add(resultPositionToOriginalLoopkup_.get(ind));
    Collections.sort(inds);
    this.currentSelectedSheet_ =  (String)selectedSheet_.getSelectedItem();
    for (int i=(inds.size()-1); i>-1; i--){
      result_.getIdentifications().get(this.currentSelectedSheet_).remove(inds.get(i).intValue());
      try {
        QuantificationResultExporter.writeResultsToExcel(selectedResultFile.getText(), result_);
        readResultFile(selectedResultFile.getText(),true);
        updateResultListSelectionTable();
      } catch (Exception e) {
        e.printStackTrace();
      }
    }
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
    this.currentSelectedSheet_ =  (String)selectedSheet_.getSelectedItem();
    return result_.getIdentifications().get(currentSelectedSheet_).get(resultPositionToOriginalLoopkup_.get(position));
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
  
  
  private class LicenseChangeListener implements ChangeListener{

//    int currentSelectedIndex_;
//    int lastSelectedIndex_;
   
    public LicenseChangeListener(){
//      currentSelectedIndex_ = 0;
    }

    public void stateChanged(ChangeEvent e)
    {
//      lastSelectedIndex_ = currentSelectedIndex_;
//      currentSelectedIndex_ = mainTabs.getSelectedIndex();
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
    lockMzRange_.setSelected(false);
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
      Hashtable<Integer,Vector<RangeColor>> rangeColors = StaticUtils.createRangeColorVOs(params,((LipidomicsTableModel)displayTable.getModel()).getMSnIdentificationName(currentSelected_),
          result_.getFaHydroxyEncoding(), result_.getLcbHydroxyEncoding(), areTheseAlex123MsnFragments());
      int threeDMsLevel = 2;
      Hashtable<Integer,Vector<String>> spectraRaw = reader_.getMsMsSpectra(params.Mz[0]-LipidomicsConstants.getMs2PrecursorTolerance(), params.Mz[0]+LipidomicsConstants.getMs2PrecursorTolerance(),-1f,-1f);

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
      if (this.show2D_.isSelected()){
        int extendedStart = (borders[0]*9)/10;
        int extendedStop = (borders[1]*21)/20;
        Lipidomics2DSpectraChromPainter spectrumPainter = new Lipidomics2DSpectraChromPainter(analyzer_,scanNrSpectrumHash, scanNrPrecursorHash, scanNrLevelHash, allRetTimes,
            peakRt, extendedStart,extendedStop,LipidomicsConstants.getMs2ChromMultiplicationFactorForInt(),LipidomicsConstants.getMs2PrecursorTolerance()*2,this,
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
          float startFloat = currentIsotopicMass-Float.parseFloat(this.displayMinusTolerance_.getText());
          float stopFloat = currentIsotopicMass+Float.parseFloat(this.displayPlusTolerance_.getText());
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
      @SuppressWarnings("unused")
      WarningMessage dlg = new WarningMessage(new JFrame(), "Error", "The MS/MS cannot be displayed: "+cgx.getMessage());
      this.displaysMs2_ = false;
      return false;
    }catch (LipidCombinameEncodingException cgx) {
      @SuppressWarnings("unused")
      WarningMessage dlg = new WarningMessage(new JFrame(), "Error", "The molecular species name cannot be decoded: "+cgx.getMessage());
      this.displaysMs2_ = false;
      return false;

    }catch (Exception cgx) {
      cgx.printStackTrace();
        @SuppressWarnings("unused")
        WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "The 3D Viewer cannot be started: "+cgx.getMessage());
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
  
  public void exportMzTab(File exportFile, short speciesType)
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
      //check if the settings type is everywhere the same, otherwise, use "relative value"
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
	        valueType = "relative value";
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
            SmallMztabMolecule molecule = MztabUtils.createSmallMztabMolecule(speciesType, summaryId, featureId, evidenceId,
	                evidenceGroupingId, maxIsotopes,analysisModule_, msRuns, originalExcelResults, molGroup, molName, resultsMol,
	                adductsSorted, expsOfGroup, faEncoding, lcbEncoding);
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
	  }
  
  /**
   * exports the mass list with the omega position information
   * @param exportFile the file to export the omega mass lists
   * @param labelInfo the information about the isotopic labels
   */
  @SuppressWarnings("unchecked")
  public void exportOmegaMasslist(File exportFile, Vector<IsotopicLabelVO> labelInfo){
    System.out.println("omega export");
    long time = System.currentTimeMillis();
    LinkedHashMap<String,Integer> classSequence = null;
    Hashtable<String,Vector<String>> correctAnalyteSequence = null;
    Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects = null;
    boolean throwError = false;
    //Hashtable<String,Hashtable<String,Hashtable<String,Vector<Double>>>> resultsStatsSection = new Hashtable<String,Hashtable<String,Hashtable<String,Vector<Double>>>>();
    if (correctOrderFile_.getText()!=null && correctOrderFile_.getText().length()>0){
      try {
        boolean ionMode = false;
        if (this.ionModeOrder_!=null && ((String)ionModeOrder_.getSelectedItem()).equalsIgnoreCase("+"))
          ionMode = true;
        @SuppressWarnings("rawtypes")
        Vector quantInfo = QuantificationThread.getCorrectAnalyteSequence(correctOrderFile_.getText(),ionMode);
        classSequence = (LinkedHashMap<String,Integer>)quantInfo.get(0);
        correctAnalyteSequence = (Hashtable<String,Vector<String>>)quantInfo.get(1);
        quantObjects = (Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>>)quantInfo.get(3);
        
//        for (String molGroup:this.heatmaps_.keySet()) {
//          HeatMapDrawing heatmap = this.heatmaps_.get(molGroup);
//          Hashtable<String,Hashtable<String,Vector<Double>>> results = heatmap.getResultValues("relative value");
//          resultsStatsSection.put(molGroup, results);
//        }

      }
      catch (Exception e) {
        throwError = true;
      }
    }else
      throwError = true;
    if (throwError) {
      new WarningMessage(new JFrame(), "Error", "For an \u03C9-RT export, you have to provide a conventional quantification file at \"\"Quant file (for correct analyte order):\"");
      return;      
    }
    try {
      Hashtable<String,Hashtable<String,Hashtable<String,Vector<AnalyteOmegaInfoVO>>>> acceptedMolecules = new Hashtable<String,Hashtable<String,Hashtable<String,Vector<AnalyteOmegaInfoVO>>>>();
      Hashtable<String,Integer> maxIsotopes = new Hashtable<String,Integer>();
      Hashtable<String,String> usedLabels = null;
      for (String molGroup:this.heatmaps_.keySet()) {
        Hashtable<String,Hashtable<String,Vector<AnalyteOmegaInfoVO>>> molsOfClass = new Hashtable<String,Hashtable<String,Vector<AnalyteOmegaInfoVO>>>();
        HeatMapDrawing heatmap = this.heatmaps_.get(molGroup);
        for (String molName: heatmap.getSelectedMoleculeNames()){
          AnalyteOmegaInfoVO omegaInfo = ifOmegaLabelExtractInfo(molGroup,molName,labelInfo);
          if (omegaInfo!=null) {
            if (!molsOfClass.containsKey(omegaInfo.getAnalyteName())) {
              molsOfClass.put(omegaInfo.getAnalyteName(), new Hashtable<String,Vector<AnalyteOmegaInfoVO>>());
            }
            usedLabels = new Hashtable<String,String>();
            for (IsotopicLabelVO label : omegaInfo.getLabels()) {
              if (usedLabels.containsKey(label.getLabelId()))
                continue;
              if (!molsOfClass.get(omegaInfo.getAnalyteName()).containsKey(label.getLabelId()))
                molsOfClass.get(omegaInfo.getAnalyteName()).put(label.getLabelId(), new Vector<AnalyteOmegaInfoVO>());
              molsOfClass.get(omegaInfo.getAnalyteName()).get(label.getLabelId()).add(omegaInfo);
              usedLabels.put(label.getLabelId(), label.getLabelId());
            }
          }
        }
        acceptedMolecules.put(molGroup, molsOfClass);
        maxIsotopes.put(molGroup, heatmap.getSelectedIsotope());
      }
      
      Hashtable<Integer,Hashtable<String,IsotopicLabelVO>> sameOmegaLabels = new  Hashtable<Integer,Hashtable<String,IsotopicLabelVO>>();
      for (IsotopicLabelVO label : labelInfo) {
        if (!sameOmegaLabels.containsKey(label.getOmegaPosition()))
          sameOmegaLabels.put(label.getOmegaPosition(), new Hashtable<String,IsotopicLabelVO>());
        sameOmegaLabels.get(label.getOmegaPosition()).put(label.getLabelId(), label);
      }

      OmegaMasslistExporter omegaExporter = new OmegaMasslistExporter(exportFile.getAbsolutePath());
      omegaExporter.export(classSequence, correctAnalyteSequence, quantObjects, acceptedMolecules, analysisModule_,sameOmegaLabels);
    }
    catch (ExportException | LipidCombinameEncodingException e) {
      e.printStackTrace();
      new WarningMessage(new JFrame(), "Error", "The \u03C9-RT mass list cannot be written: "+e.getMessage());
    }
    long usedTime = (System.currentTimeMillis()-time)/1000l;
    System.out.println("Used time: "+(usedTime/60l)+" minutes "+usedTime%60l+" seconds");
  }
  
  
  /**
   * extracts information on whether a detection is an isotopically labeled species and which labels have been applied
   * @param className name of the lipid class
   * @param name the name of the molecule
   * @param labelInfo the value objects containg information about a label
   * @return an AnalyteOmegaInfoVO if the detection is an isotopically labeled species, null otherwise
   */
  private AnalyteOmegaInfoVO ifOmegaLabelExtractInfo(String className, String name, Vector<IsotopicLabelVO> labelInfo) {
    boolean isLabel = false;
    //if it does not start with a label indicator, we can immediately return false;
    for (IsotopicLabelVO info : labelInfo) {
      if (name.startsWith(info.getLabelId())) {
        isLabel = true;
        break;
      }
    }
    if (!isLabel)
      return null; 
    //the avoid false positive label identifications by molecules whose name simply starts with the same letter as a label, we have to make a further check
    String analyze = new String(name);
    boolean isALabel = true;
    boolean foundLabel;
    Vector<IsotopicLabelVO> detectedLabels = new Vector<IsotopicLabelVO>();
    while (isALabel) {
      foundLabel = false;
      for (IsotopicLabelVO info : labelInfo) {
        if (analyze.startsWith(info.getLabelId())) {
          foundLabel = true;
          analyze = analyze.substring(info.getLabelId().length());
          detectedLabels.add(info);
          break;
        }
      }
      if (!foundLabel)
        isALabel = false;
    }
    //TODO: this does not work if no RT-grouping is set
    String analyteName = analyze.substring(0,analyze.lastIndexOf("_"));
    //this is for sphingolipids
    if (analyze.length()>0 && analysisModule_.getLcbHydroxyEncoding()!=null && analysisModule_.getLcbHydroxyEncoding().containsValue(analyze.substring(0,1)))
      analyze = analyze.substring(1);
    if (analyze.length()==0)
      return null;
    if (Character.isDigit(analyze.toCharArray()[0])) {
      boolean isMSnEvidence = false;
      ResultAreaVO resultVO = null;
      for (String exp : this.analysisModule_.getExpNamesInSequence()) {
        resultVO = analysisModule_.getResultAreaVO(className,name,exp);
        if (resultVO!=null && resultVO.isMsnEvidenceThere()) {
          isMSnEvidence = true;
          break;
        }
      }
      if (!isMSnEvidence)
        return null;
      boolean oneDoubleBond = false;
      if (analyteName.indexOf(":")!=-1 && analyteName.substring(analyteName.indexOf(":")+1).equalsIgnoreCase("1"))
        oneDoubleBond = true;
      return new AnalyteOmegaInfoVO(analyteName,name,detectedLabels,oneDoubleBond);
    }else
      return null;
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
  
  
  private void cleanupResultView(){
    if (analysisModule_!=null){
      analysisModule_.cleanup();
      if (groupDisplayNamesLookup_!=null)groupDisplayNamesLookup_.clear();
      groupDisplayNamesLookup_ = null;
      if (expDisplayNamesLookup_!=null)expDisplayNamesLookup_.clear();
      expDisplayNamesLookup_ = null;
      if (heatmaps_!=null){
        for (HeatMapDrawing map : heatmaps_.values())map.cleanup();
        heatmaps_.clear();
        heatmaps_ = null;
      }
      if (groupHeatmaps_!=null){
        for (HeatMapDrawing map : groupHeatmaps_.values())map.cleanup();
        groupHeatmaps_.clear();
        groupHeatmaps_ = null;
      }
      System.gc();
    }
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
    Hashtable<Integer,Vector<RangeColor>> rangeColors = StaticUtils.createRangeColorVOs(param, ((LipidomicsTableModel)displayTable.getModel()).getMSnIdentificationName(currentSelected_),
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
    lockMzRange_.setSelected(false);
    if(newPrec ==  true){
      //PAINTS THE NEW SPEKTRA WITH THE LPS
      float m2dGain = spectrumPainter_.getM2dGain();   
      spectrumPanel_.remove(spectrumPainter_);   
      //int msLevel = 2;
      Hashtable<Integer,Vector<String>> spectraRaw = reader_.getMsMsSpectra(params.Mz[0]-LipidomicsConstants.getMs2PrecursorTolerance(), params.Mz[0]+LipidomicsConstants.getMs2PrecursorTolerance(), -1f, -1f);
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
          peakRt, extendedStart,extendedStop,LipidomicsConstants.getMs2ChromMultiplicationFactorForInt(),LipidomicsConstants.getMs2PrecursorTolerance()*2,this,
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
      float startFloat = currentIsotopicMass-Float.parseFloat(this.displayMinusTolerance_.getText());
      float stopFloat = currentIsotopicMass+Float.parseFloat(this.displayPlusTolerance_.getText());
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
      if (!this.lockMzRange_.isSelected() || this.viewer_==null || lockRangeUpdateRequired_){
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
    if (params!=null && this.show2D_.isSelected()) {
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
    boolean selected = lockMzRange_.isSelected();
    displayMinusTolerance_.setEnabled(!selected);
    displayPlusTolerance_.setEnabled(!selected);
    displayMzStart_.setEnabled(selected);
    displayMzStop_.setEnabled(selected);
    if (lockMzRange_.isSelected() && params_!=null && (displayMzStart_.getText()==null || displayMzStart_.getText().length()==0)) {
      float currentIsotopicMass = params_.Mz[0];
      displayMzStart_.setText(Calculator.FormatNumberToString((double)(currentIsotopicMass-Float.parseFloat(this.displayMinusTolerance_.getText())),2d));
      displayMzStop_.setText(Calculator.FormatNumberToString((double)(currentIsotopicMass+Float.parseFloat(this.displayPlusTolerance_.getText())),2d));
    }
    if (lockMzRange_.isSelected())
      lockRangeUpdateRequired_ = true;
  }
}
  
 