package at.tugraz.genome.lda.analysis;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.GradientPaint;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.image.BufferedImage;
import java.util.Hashtable;
import java.util.Vector;

import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.maspectras.graphics.MacOSResizer;

/**
 * 
 * @author Jürgen Hartler
 * @author Leonida M. Lamp
 *
 */
public abstract class HeatMapGenerator {

  public final static int NUMBER_OF_COLOR_GROUPS = 9;
  Color borderColor_=new Color(0,0,0);
  
  protected transient BufferedImage expressionImage;
  protected Vector<String> sampleNames_;
  protected Vector<String> analyteNames_;
  
  protected int defaultNrOfGradientPixels_ = 10000;

  protected final static Color DEFAULT_COLOR0 = new Color(000,040,180);
  protected final static Color DEFAULT_COLOR1 = new Color(000,153,102);
  protected final static Color DEFAULT_COLOR2 = new Color(255,255,000);
  protected final static Color DEFAULT_COLOR3 = new Color(255,255,000);
  protected final static Color DEFAULT_COLOR4 = new Color(255,000,000);
  private final static String LABEL_DB_CONFLICT = "\u03C9-DB assignment conflict";
  private final static String LABEL_DB_ASSIGNED = "\u03C9-DB assigned";
  
  protected BufferedImage gradient_;
  
  protected int heatRectWidth_ = 20;
  protected int heatRectHeight_ = 10;
  protected int borderCorrection_ = 1;
  
  protected int annotationSpace_ = 200;
  protected int sampleSpace_ = 80;
  protected int pictureIndent_ = 5;
  protected int gradientHeight_ = 10;
  protected int analyteNameSpace_ = 3;
//  protected int doubleBondNameSpace_ = 5;
  
  protected int analyteTextWidthMax_ = 0;
  protected int sampleNameStart_ = 0;
  
  protected int doubleBondTextWidthMax_ = 0;
  
  protected Hashtable<String,Hashtable<String,Color>> attentionValues_;
  
  int partCorrection_ = 0;
  
  float[][] data_;
  
  protected abstract BufferedImage createThreeColGradientImage(int sizeX);
  
  protected abstract void paintHeader(Graphics2D expressionGraphics,int imageSizeX);
  
  protected abstract int getColorForValue(float value);
  
  protected abstract String getSampleNameToDisplay(String sampleName);
  
  protected abstract boolean isOSMacAndJavaLookAndFeel();
  
  protected abstract Color getLegendColor();
  
  protected abstract Color getBackgroundColor();
  
  
//  protected void setInputValues(float[][] data, Vector<String> sampleNames, Vector<String> analyteNames, Hashtable<String,Hashtable<String,Color>> attentionValues){
//    this.data_ = data;
//    this.sampleNames_ = sampleNames;
//    this.analyteNames_ = analyteNames;  
//    this.attentionValues_ = attentionValues;
//  }
//  
//  public void init() {
//    this.initGradient();
//  }
//  
//  private void initGradient(){
//    this.calculateMinMax();
//    //gradient_ = HeatMapGenerator.createDefaultGradientImage(this.defaultNrOfGradientPixels_);
//    gradient_ = this.createThreeColGradientImage(this.defaultNrOfGradientPixels_);
//    
//  }
  
  public BufferedImage createImage() {
    return this.createImage(null);
  }
  
  public BufferedImage createImage(Graphics2D extGraphics) {
    BufferedImage dummyImage = new BufferedImage(1000,1000,BufferedImage.TYPE_3BYTE_BGR);
    Font descriptionFont = new Font("Dialog",Font.PLAIN, 9);
    Graphics2D dummyGraphics = (Graphics2D)dummyImage.getGraphics();
    dummyGraphics.setFont(descriptionFont);
    FontMetrics descriptionFontMetrics = dummyGraphics.getFontMetrics();
    int textWidth = 0;
    for (String sampleName : sampleNames_){
      String  toDisplay = getSampleNameToDisplay(sampleName);
      int width = descriptionFontMetrics.stringWidth(toDisplay);
      if (width>textWidth){
        textWidth = width;
      }  
    }
    int annotationWidth = 0;
    for (String analyteName : analyteNames_){
      int width = descriptionFontMetrics.stringWidth(analyteName);
      if (width>annotationWidth){
        annotationWidth = width;
      }  
    }
//    int doubleBondTextWidthMax = 0;
//    if (LipidParameterSet.isOmegaInformationAvailable()) {
//      doubleBondTextWidthMax = descriptionFontMetrics.stringWidth(LABEL_DB_CONFLICT);
//    }
    try {
      this.sampleNameStart_ = this.pictureIndent_+this.gradientHeight_+descriptionFontMetrics.getHeight();
      this.sampleSpace_ = textWidth+sampleNameStart_;
      this.annotationSpace_ = annotationWidth+2*this.pictureIndent_;
//      this.annotationSpace_ = annotationWidth+doubleBondTextWidthMax+2*this.pictureIndent_;
      Graphics2D expressionGraphics;
      int xcorrection=0;
      int ycorrection=0;
      int xSpaceForMap = this.heatRectWidth_*this.sampleNames_.size()+borderCorrection_;
      int imageSizeX = xSpaceForMap+/*(this.getNumberOfParts()-1)*10*/+xcorrection+annotationSpace_+2*pictureIndent_;
      expressionImage=new BufferedImage(imageSizeX,
                                         this.heatRectHeight_*(this.analyteNames_.size())+borderCorrection_+ycorrection+sampleSpace_+2*pictureIndent_,
                                         BufferedImage.TYPE_3BYTE_BGR);
       expressionGraphics=(Graphics2D)expressionImage.getGraphics();
       if (extGraphics!=null)
         expressionGraphics=extGraphics;
       expressionGraphics.setColor(getBackgroundColor());
       expressionGraphics.fillRect(0, 0, expressionImage.getWidth(), expressionImage.getHeight());
       this.paintHeader(expressionGraphics, imageSizeX-annotationSpace_);
       this.paintSampleNames(expressionGraphics,this.pictureIndent_ /*imageSizeX-pictureIndent_-xSpaceForMap-1*/,sampleSpace_);
       paintExpressionImage(expressionGraphics, this.pictureIndent_/*imageSizeX-pictureIndent_-xSpaceForMap-1*/, sampleSpace_/*, geneOrder, sampleOrder*/);

       this.paintAnalyteNames(expressionGraphics,this.pictureIndent_+xSpaceForMap,sampleSpace_);
//       if (LipidParameterSet.isOmegaInformationAvailable()) {
//         this.paintLipidDoubleBondNames(expressionGraphics,this.pictureIndent_+xSpaceForMap+analyteTextWidthMax_+doubleBondNameSpace_,sampleSpace_);
//       }
       
       } catch (OutOfMemoryError e) {
       }
       return expressionImage;      
 }
  

  private void paintExpressionImage(Graphics2D expressionGraphics, int xoffset,
      int yoffset)
  {
    float dummyfloat;
    expressionGraphics.setColor(borderColor_);
//    expressionGraphics.fillRect(xoffset, yoffset, expressionImage.getWidth(),
//        expressionImage.getHeight());
    

    for (int i = 0; i < this.analyteNames_.size(); i++) {
      partCorrection_ = 0;
      for (int j = 0; j < this.sampleNames_.size(); j++) {
        dummyfloat = this.getValue(j, i);
        Color dummycolor = new Color(getColorForValue(dummyfloat));
        expressionGraphics.setColor(dummycolor);
        expressionGraphics.fillRect(this.heatRectWidth_ * j + 1
            + xoffset + partCorrection_, this.heatRectHeight_ * i + 1
            + yoffset, this.heatRectWidth_ - 1, this.heatRectHeight_ - 1);
        if (attentionValues_.containsKey(sampleNames_.get(j)) && attentionValues_.get(sampleNames_.get(j)).containsKey(analyteNames_.get(i))){
          paintAttentionRectangle(expressionGraphics,attentionValues_.get(sampleNames_.get(j)).get(analyteNames_.get(i)),i,j,xoffset,yoffset);
        }              
      }
    }
  }

  
  /**
   * this method paints an unfilled rectangle around the heat map rectangle containing the expression color
   * @param g2 Graphics2 object for painting
   * @param color the color to be used for painting the unfilled rectangle
   * @param row the row number in the expression image
   * @param column the column number in the expression image
   * @param xoffset the offset in x direction of the heat map within the image
   * @param yoffset the offset in y direction of the heat map within the image
   */
  public void paintAttentionRectangle(Graphics2D g2, Color color, int row, int column, int xoffset, int yoffset){
    g2.setColor(attentionValues_.get(sampleNames_.get(column)).get(analyteNames_.get(row)));
    g2.drawRect(this.heatRectWidth_ * column + 1
        + xoffset + partCorrection_, this.heatRectHeight_ * row + 1
        + yoffset, this.heatRectWidth_ - 2, this.heatRectHeight_ - 2);
  }

  
  private void paintSampleNames(Graphics2D expressionGraphics, int xoffset, int sampleSpace_){
    Font bigFont = new Font("SansSerif", Font.PLAIN, 200);
    expressionGraphics.setFont(bigFont);
    FontMetrics fm2 = expressionGraphics.getFontMetrics();

    Font descriptionFont = new Font("Dialog",Font.PLAIN, 9);
    expressionGraphics.setColor(getLegendColor());
    
    expressionGraphics.setFont(descriptionFont);
    FontMetrics descriptionFontMetrics = expressionGraphics.getFontMetrics();
    int textHeight = descriptionFontMetrics.getHeight();
    double rotationAngle = -Math.PI/2.0;
    if (!isOSMacAndJavaLookAndFeel())
      expressionGraphics.rotate(rotationAngle);
    for (int i=0;i!=this.sampleNames_.size();i++){
      String sampleName = sampleNames_.get(i);
      String toDisplay = this.getSampleNameToDisplay(sampleName);
      if (isOSMacAndJavaLookAndFeel()){       
        BufferedImage image = new BufferedImage(fm2.stringWidth(toDisplay), fm2.getHeight(),
            BufferedImage.TYPE_INT_ARGB);
        Graphics2D g22 = (Graphics2D) image.getGraphics();
        g22.setColor(getLegendColor());
        g22.setFont(bigFont);
        g22.drawString(toDisplay,0,fm2.getAscent());
        int sw2 = fm2.stringWidth(toDisplay);
        int swXExtension2 = (int)(((double)sw2)*Math.cos(-rotationAngle));
        int swYExtension2 = (int)(((double)sw2)*Math.sin(-rotationAngle));
        int shYExtension2 = (int)(((double)fm2.getHeight())*Math.cos(-rotationAngle));
        int shXExtension2 = (int)(((double)fm2.getHeight())*Math.sin(-rotationAngle));
      
        BufferedImage image2 = new BufferedImage(swXExtension2+shXExtension2, swYExtension2+shYExtension2,
            BufferedImage.TYPE_INT_ARGB);
        Graphics2D g23 = (Graphics2D) image2.getGraphics();
        g23.rotate(rotationAngle);
        int normalStringXCoordinate2 = 0;
        int normalStringYCoordinate2 = swYExtension2;
        double coordinateRadius2 = Math.sqrt(Math.pow(0, 2)+Math.pow(normalStringYCoordinate2, 2));
        double angle2 = Math.atan((double)normalStringYCoordinate2/(double)normalStringXCoordinate2);
        g23.drawImage(image,(int)(coordinateRadius2*Math.cos(angle2-rotationAngle)),(int)(coordinateRadius2*Math.sin(angle2-rotationAngle)),null);
        g23.rotate(-rotationAngle);   
        //expressionGraphics.drawImage(image2.getScaledInstance(descriptionFontMetrics.getHeight(), descriptionFontMetrics.stringWidth(toDisplay), 0),this.heatRectWidth_ * i +this.heatRectWidth_/2,sampleSpace_-2-descriptionFontMetrics.stringWidth(toDisplay),null);
        image2 = MacOSResizer.resizeTrick(image2, descriptionFontMetrics.getHeight(), descriptionFontMetrics.stringWidth(toDisplay));
        expressionGraphics.drawImage(image2,this.heatRectWidth_ * i +this.heatRectWidth_/2,sampleSpace_-2-descriptionFontMetrics.stringWidth(toDisplay),null);
      }else{
        expressionGraphics.drawString(toDisplay,-sampleSpace_+2,this.heatRectWidth_ * i +this.heatRectWidth_/2+textHeight/2
            + xoffset);
      }
    }
    if (!isOSMacAndJavaLookAndFeel())
      expressionGraphics.rotate(-rotationAngle);

  }
  
  private void paintAnalyteNames(Graphics2D expressionGraphics, int xoffset, int yoffset){
    Font descriptionFont = new Font("Dialog",Font.PLAIN, 9);
    expressionGraphics.setColor(getLegendColor());
    expressionGraphics.setFont(descriptionFont);
    FontMetrics descriptionFontMetrics = expressionGraphics.getFontMetrics();
    int textHeight = descriptionFontMetrics.getHeight();
    for (int i=0;i!=this.analyteNames_.size();i++){
      String analyteName = this.analyteNames_.get(i);
      int textWidth = descriptionFontMetrics.stringWidth(analyteName);
      if (textWidth>analyteTextWidthMax_)
        analyteTextWidthMax_ = textWidth;
      expressionGraphics.drawString(analyteName,xoffset+analyteNameSpace_,(this.heatRectHeight_*i) + 2 + yoffset+textHeight/2);
    }
  }
  
//  private void paintLipidDoubleBondNames(Graphics2D expressionGraphics, int xoffset, int yoffset){
//    Font descriptionFont = new Font("Dialog",Font.PLAIN, 9);
//    expressionGraphics.setColor(getLegendColor());
//    expressionGraphics.setFont(descriptionFont);
//    FontMetrics descriptionFontMetrics = expressionGraphics.getFontMetrics();
//    doubleBondTextWidthMax_ = descriptionFontMetrics.stringWidth(LABEL_DB_CONFLICT);
//    int textHeight = descriptionFontMetrics.getHeight();
//    for (int i=0;i<this.analyteNames_.size();i++){
//      String analyteName = this.analyteNames_.get(i);
//      String doubleBondPositionInformation = null;
//      int doubleBondPositionAvailability = this.doubleBondPositionAvailability_.get(analyteName);
//      switch (doubleBondPositionAvailability) {
//        case 0:
//          doubleBondPositionInformation = "";
//          break;
//        case 1:
//          doubleBondPositionInformation = LABEL_DB_CONFLICT;
//          break;
//        case 2:
//          doubleBondPositionInformation = LABEL_DB_ASSIGNED;
//          break;
//        default:
//          break;
//      }
//      expressionGraphics.drawString(doubleBondPositionInformation,xoffset+doubleBondNameSpace_,(this.heatRectHeight_*i) + 2 + yoffset+textHeight/2);
//    }
//  }
  
  protected float[] mapColorToScale(int rgb1){
    float[] color = new float[3];
    color[0] = (float)((rgb1 >> 16) & 0xFF)/255f;
    color[1] = (float)((rgb1 >>  8) & 0xFF)/255f;
    color[2] = (float)((rgb1 >>  0) & 0xFF)/255f;
    return color;    
  }
  
  
  
  protected static BufferedImage createDefaultGradientImage(int sizeX) {
    BufferedImage gradient = new BufferedImage(sizeX, 1, BufferedImage.TYPE_3BYTE_BGR);
    Graphics2D graphics = gradient.createGraphics();
//    GradientPaint gp1 = new GradientPaint(0, 0, HeatMapGenerator.DEFAULT_COLOR1, (sizeX-1)/15, 0,HeatMapGenerator.DEFAULT_COLOR2);
//    GradientPaint gp2 = new GradientPaint((sizeX-1)/15, 0, HeatMapGenerator.DEFAULT_COLOR3,  sizeX-1, 0, HeatMapGenerator.DEFAULT_COLOR4);
    GradientPaint gp1 = new GradientPaint(0, 0, HeatMapGenerator.DEFAULT_COLOR1, (sizeX-1)/2, 0,HeatMapGenerator.DEFAULT_COLOR2);
    GradientPaint gp2 = new GradientPaint((sizeX-1)/2, 0, HeatMapGenerator.DEFAULT_COLOR3,  sizeX-1, 0, HeatMapGenerator.DEFAULT_COLOR4);
    graphics.setPaint(gp1);
    graphics.drawRect(0, 0, (sizeX-1)/2, 1);
    graphics.setPaint(gp2);
    graphics.drawRect((sizeX-1)/2, 0, sizeX-1, 1);
    graphics.dispose();    
    return gradient;
  }
  
  public float getValue(int x, int y) {
    return data_[x][y];
 }
  
  public int getExpressionImageXStart(){
    return this.pictureIndent_;
  }

  public int getExpressionImageXEnd(){
    return this.pictureIndent_+this.heatRectWidth_*sampleNames_.size();
  }
  
  public int getExpressionImageYStart(){
    return sampleSpace_;
  }

  public int getExpressionImageYEnd(){
    return this.sampleSpace_+this.heatRectHeight_*this.analyteNames_.size();
  }
  
  public Vector<String> getCellName(int x, int y){
    Vector<String> cellName = new Vector<String>();
    String sampleName = null;
    String analyteName = null;
    int[] cellPostion = this.getCellPosition(x, y);
    if (cellPostion[0]>=0&&cellPostion[1]>=0){
      sampleName = sampleNames_.get(cellPostion[0]);
      analyteName = analyteNames_.get(cellPostion[1]);
    }
    cellName.add(sampleName);
    cellName.add(analyteName);
    return cellName;
  }
  
  public int[] getCellPosition(int x, int y){
    int[] cellPosition = new int[2];
    cellPosition[0] = -1;
    cellPosition[1] = -1;
    if (this.getExpressionImageXStart()<=x&&x<this.getExpressionImageXEnd()&&
        this.getExpressionImageYStart()<=y&&y<this.getExpressionImageYEnd()){
      cellPosition[0] =this.calculateXPosition(x);
      cellPosition[1] = this.calculateYPosition(y);;
    }
    return cellPosition;
  }
  
  public Rectangle getRectangleForCell(int x, int y){
    Rectangle rect = null;
    if (this.getExpressionImageXStart()<=x&&x<this.getExpressionImageXEnd()&&
        this.getExpressionImageYStart()<=y&&y<this.getExpressionImageYEnd()){
      int cellColumn = this.calculateXPosition(x);
      int cellRow = this.calculateYPosition(y);
      rect = getRectangleForCellByRowAndColumn(cellRow, cellColumn);
    }
    return rect;
  }
  
  /**
   * returns a rectangle of a heat map cell defined by row and column number
   * @param cellRow the row number of the cell
   * @param cellColumn the column number of the cell
   * @return rectangle of a heat map cell defined by row and column number
   */
  public Rectangle getRectangleForCellByRowAndColumn(int cellRow, int cellColumn){
    return new Rectangle(this.heatRectWidth_ * cellColumn+getExpressionImageXStart(),
        this.heatRectHeight_ * cellRow+getExpressionImageYStart(),this.heatRectWidth_,
        this.heatRectHeight_);
  }
  
  public int getRowNameStart(){
    return this.getExpressionImageXEnd()+this.analyteNameSpace_;
  }
  
  public int getRowNameEnd(){
    return this.getRowNameStart()+this.analyteTextWidthMax_;
  }
  
  public String getRowName(int x, int y){
    String analyteName = null;
    if (this.getRowNameStart()<=x&&x<this.getRowNameEnd()&&
        this.getExpressionImageYStart()<=y&&y<this.getExpressionImageYEnd()){
    	analyteName = analyteNames_.get(this.calculateYPosition(y));
    }
    return analyteName;
  }
  
  public Rectangle getRectangleForRowName(String rowName){
    int rowPosition = -1;
    int count = 0;
    for (String name : this.analyteNames_){
      if (rowName.equalsIgnoreCase(name)){
        rowPosition = count;
        break;
      }
      count++;
    }
    Rectangle rect = null;
    if (rowPosition>=0){
      rect = new Rectangle(this.getRowNameStart(),this.heatRectHeight_ * rowPosition+getExpressionImageYStart(),
          this.analyteTextWidthMax_+this.analyteNameSpace_, this.heatRectHeight_);
    }
    return rect;
  }
  
//  public int getLipidDoubleBondRowNameStart(){
//    return this.getExpressionImageXEnd()+this.analyteNameSpace_+analyteTextWidthMax_+doubleBondNameSpace_;
//  }
//  
//  public int getLipidDoubleBondRowNameEnd(){
//    return this.getLipidDoubleBondRowNameStart()+this.doubleBondTextWidthMax_;
//  }
  
//  public String[] getLipidDoubleBondRowName(int x, int y){
//    String[] doubleBondInfo = new String[2];
//    if (this.getLipidDoubleBondRowNameStart()<=x&&x<this.getLipidDoubleBondRowNameEnd()&&
//        this.getExpressionImageYStart()<=y&&y<this.getExpressionImageYEnd()){
//      doubleBondInfo[0] = analyteNames_.get(this.calculateYPosition(y));
//      int availability = doubleBondPositionAvailability_.get(doubleBondInfo[0]);
//      switch (availability) {
//        case 0:
//          doubleBondInfo[1] = "";
//          break;
//        case 1:
//          doubleBondInfo[1] = LABEL_DB_CONFLICT;
//          break;
//        case 2:
//          doubleBondInfo[1] = LABEL_DB_ASSIGNED;
//          break;
//        default:
//          break;
//      }
//    }
//    return doubleBondInfo;
//  }
  
//  public Rectangle getRectangleForLipidDoubleBondRowName(String rowName){
//    int rowPosition = -1;
//    int count = 0;
//    for (String name : this.analyteNames_){
//      if (rowName.equalsIgnoreCase(name)){
//        rowPosition = count;
//        break;
//      }
//      count++;
//    }
//    Rectangle rect = null;
//    if (rowPosition>=0){
//      rect = new Rectangle(this.getLipidDoubleBondRowNameStart(),this.heatRectHeight_ * rowPosition+getExpressionImageYStart(),
//          this.doubleBondTextWidthMax_+this.doubleBondNameSpace_, this.heatRectHeight_);
//    }
//    return rect;
//  }
  
  private int calculateXPosition(int x){
    return (x-this.getExpressionImageXStart())/heatRectWidth_;
  }  
  
  private int calculateYPosition(int y){
    return (y-this.getExpressionImageYStart())/heatRectHeight_;
  }

  public int getColumnNameStart(){
    return this.sampleNameStart_-2;
  }
  
  public int getColumnNameEnd(){
    return this.sampleSpace_-2;
  }

  public String getColumnName(int x, int y){
    String sampleName = null;
    if (this.getExpressionImageXStart()<=x&&x<this.getExpressionImageXEnd()&&
        this.getColumnNameStart()<=y&&y<this.getColumnNameEnd()){
      sampleName = sampleNames_.get(this.calculateXPosition(x));
    }
    return sampleName;
  }
  
  public Rectangle getRectangleForColumnName(String columnName){
    int columnPosition = -1;
    int count = 0;
    for (String name : this.sampleNames_){
      if (columnName.equalsIgnoreCase(name)){
        columnPosition = count;
        break;
      }
      count++;
    }
    Rectangle rect = null;
    if (columnPosition>=0){
      rect = new Rectangle(this.heatRectWidth_ * columnPosition+getExpressionImageXStart(),this.getColumnNameStart(),
          this.heatRectWidth_, this.getColumnNameEnd()-this.getColumnNameStart());
    }
    return rect;
  }
  
  /**
   * returns the color of a heat map cell defined by row and column number
   * @param cellRow the row number of the cell
   * @param cellColumn the column number of the cell
   * @return color of a heat map cell defined by row and column number
   */
  public Color getColorForCell(int cellRow, int cellColumn){
    return new Color(getColorForValue(this.getValue(cellColumn, cellRow)));
  }
  
}
