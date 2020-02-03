package at.tugraz.genome.vos;

import java.util.LinkedHashMap;

public class Mix2Standards
{
  
  
  public static LinkedHashMap<String,ReferenceInfoVO> getSphBaseStandards(){
    LinkedHashMap<String,ReferenceInfoVO> standards = new LinkedHashMap<String,ReferenceInfoVO>();
    // the first boolean is if the standard shall be used in the evaluation, and the second if a position is detectable
    standards.put("d14:1", new ReferenceInfoVO("d14:1",new double[]{1.2d},true,false));
    standards.get("d14:1").addAdduct("H", 244.227105577d);
    standards.get("d14:1").addAdduct("-OH", 226.216540877d);
    standards.get("d14:1").addAdduct("Na", 266.209051277d);
    standards.put("m18:0", new ReferenceInfoVO("m18:0",new double[]{4.2d},false,false));
    standards.get("m18:0").addAdduct("H", 286.310441297d);
    standards.get("m18:0").addAdduct("-OH", 268.299876597d);
    standards.get("m18:0").addAdduct("Na", 308.292386997d);
    standards.put("d18:0", new ReferenceInfoVO("d18:0",new double[]{3.4d},true,false));
    standards.get("d18:0").addAdduct("H", 302.305355927d);
    standards.get("d18:0").addAdduct("-OH", 284.294791227d);
    standards.get("d18:0").addAdduct("Na", 324.287301627d);
    standards.put("d18:2", new ReferenceInfoVO("d18:2",new double[]{2.3d},true,false));
    standards.get("d18:2").addAdduct("H", 298.274055787d);
    standards.get("d18:2").addAdduct("-OH", 280.263491087d);
    standards.get("d18:2").addAdduct("Na", 320.256001487d);
    standards.put("d20:0", new ReferenceInfoVO("d20:0",new double[]{5.7d},true,false));
    standards.get("d20:0").addAdduct("H", 330.336656067d);
    standards.get("d20:0").addAdduct("-OH", 312.326091367d);
    standards.get("d20:0").addAdduct("Na", 352.318601767d);
    standards.put("d20:1", new ReferenceInfoVO("d20:1",new double[]{5.1d},true,false));
    standards.get("d20:1").addAdduct("H", 328.321005997d);
    standards.get("d20:1").addAdduct("-OH", 310.310441297d);
    standards.get("d20:1").addAdduct("Na", 350.302951697d);    
    return standards;
  }

  public static LinkedHashMap<String,ReferenceInfoVO> getCerStandards(){
    LinkedHashMap<String,ReferenceInfoVO> standards = new LinkedHashMap<String,ReferenceInfoVO>();
    // the first boolean is if the standard shall be used in the evaluation, and the second if a position is detectable
    standards.put("d24:0", new ReferenceInfoVO("d18:0/n6:0",new double[]{9.6d},true,false));
    standards.get("d24:0").addAdduct("-H", 398.363968307d);
    standards.get("d24:0").addAdduct("HCOO", 444.369447637d);
    standards.get("d24:0").addAdduct("Cl", 434.340646607d);
    standards.get("d24:0").addAdduct("H", 400.378520907d);
    standards.get("d24:0").addAdduct("-OH", 382.367956207d);
    standards.get("d24:0").addAdduct("Na", 422.360466607d);
    standards.put("t30:0", new ReferenceInfoVO("d18:0/h12:0",new double[]{16.3d,16.9d},true,false));
    standards.get("t30:0").addAdduct("-H", 498.452783357d);
    standards.get("t30:0").addAdduct("HCOO", 544.458262687d);
    standards.get("t30:0").addAdduct("Cl", 534.429461657d);
    standards.get("t30:0").addAdduct("H", 500.467335957d);
    standards.get("t30:0").addAdduct("-OH", 482.456771257d);
    standards.get("t30:0").addAdduct("Na", 522.449281657d);
    standards.put("m30:1", new ReferenceInfoVO("m18:1/n12:0",new double[]{19.3d},false,false));
    standards.get("m30:1").addAdduct("-H", 464.447304027d);
    standards.get("m30:1").addAdduct("HCOO", 510.452783357d);
    standards.get("m30:1").addAdduct("Cl", 500.423982327d);
    standards.get("m30:1").addAdduct("H", 466.461856627d);
    standards.get("m30:1").addAdduct("-OH", 448.451291927d);
    standards.get("m30:1").addAdduct("Na", 488.443802327d);
    standards.put("d34:1", new ReferenceInfoVO("d18:1/n16:0",new double[]{22.1d},true,false));
    standards.get("d34:1").addAdduct("-H", 536.504818937d);
    standards.get("d34:1").addAdduct("HCOO", 582.510298267d);
    standards.get("d34:1").addAdduct("Cl", 572.4814972377d);
    standards.get("d34:1").addAdduct("H", 538.519371537d);
    standards.get("d34:1").addAdduct("-OH", 520.508806837d);
    standards.get("d34:1").addAdduct("Na", 560.501317237d);  
    standards.put("d35:1", new ReferenceInfoVO("d18:1/n17:0",new double[]{23.2d},true,false));
    standards.get("d35:1").addAdduct("-H", 550.520469007d);
    standards.get("d35:1").addAdduct("HCOO", 596.525948337d);
    standards.get("d35:1").addAdduct("Cl", 586.497147307d);
    standards.get("d35:1").addAdduct("H", 552.535021607d);
    standards.get("d35:1").addAdduct("-OH", 534.524456907d);
    standards.get("d35:1").addAdduct("Na", 574.516967307d);
    standards.put("d36:1", new ReferenceInfoVO("d18:1/n18:0",new double[]{24.3d},true,false));
    standards.get("d36:1").addAdduct("-H", 564.536119077d);
    standards.get("d36:1").addAdduct("HCOO", 610.541598407d);
    standards.get("d36:1").addAdduct("Cl", 600.512797377d);
    standards.get("d36:1").addAdduct("H", 566.550671677d);
    standards.get("d36:1").addAdduct("-OH", 548.540106977d);
    standards.get("d36:1").addAdduct("Na", 588.532617377d);
    standards.put("t36:2", new ReferenceInfoVO("d18:1/h18:1",new double[]{21.6d},true,false));
    standards.get("t36:2").addAdduct("-H", 578.515383637d);
    standards.get("t36:2").addAdduct("HCOO", 624.520862967d);
    standards.get("t36:2").addAdduct("Cl", 614.492061937d);
    standards.get("t36:2").addAdduct("H", 580.529936237d);
    standards.get("t36:2").addAdduct("-OH", 562.519371537d);
    standards.get("t36:2").addAdduct("Na", 602.511881937d);
    standards.put("d38:1", new ReferenceInfoVO("d18:1/n20:0",new double[]{26.3d},true,false));
    standards.get("d38:1").addAdduct("-H", 592.567419217d);
    standards.get("d38:1").addAdduct("HCOO", 638.572898547d);
    standards.get("d38:1").addAdduct("Cl", 628.544097517d);
    standards.get("d38:1").addAdduct("H", 594.581971817d);
    standards.get("d38:1").addAdduct("-OH", 576.571407117d);
    standards.get("d38:1").addAdduct("Na", 616.563917517d);
    standards.put("d42:1", new ReferenceInfoVO("d18:0/n24:1",new double[]{28.6d},true,false));
    standards.get("d42:1").addAdduct("-H", 648.630019497d);
    standards.get("d42:1").addAdduct("HCOO", 694.635498827d);
    standards.get("d42:1").addAdduct("Cl", 684.606697797d);
    standards.get("d42:1").addAdduct("H", 650.644572097d);
    standards.get("d42:1").addAdduct("-OH", 632.634007397d);
    standards.get("d42:1").addAdduct("Na", 672.626517797d);
    standards.put("t42:1", new ReferenceInfoVO("d18:1/h24:0",new double[]{28.7d},true,false));
    standards.get("t42:1").addAdduct("-H", 664.624934127d);
    standards.get("t42:1").addAdduct("HCOO", 710.630413457d);
    standards.get("t42:1").addAdduct("Cl", 700.601612427d);
    standards.get("t42:1").addAdduct("H", 666.639486727d);
    standards.get("t42:1").addAdduct("-OH", 648.628922027d);
    standards.get("t42:1").addAdduct("Na", 688.621432427d);
    standards.put("t48:1", new ReferenceInfoVO("d18:1/h30:0",new double[]{29.6d},true,false));
    standards.get("t48:1").addAdduct("-H", 748.718834547d);
    standards.get("t48:1").addAdduct("HCOO", 794.724313877d);
    standards.get("t48:1").addAdduct("Cl", 784.6955128477d);
    standards.get("t48:1").addAdduct("H", 750.733387147d);
    standards.get("t48:1").addAdduct("-OH", 732.722822447d);
    standards.get("t48:1").addAdduct("Na", 772.715332847d);
    return standards;
  }
  
  
  public static LinkedHashMap<String,ReferenceInfoVO> getCer1PStandards(){
    LinkedHashMap<String,ReferenceInfoVO> standards = new LinkedHashMap<String,ReferenceInfoVO>();
    // the first boolean is if the standard shall be used in the evaluation, and the second if a position is detectable
    standards.put("d25:1", new ReferenceInfoVO("d17:1/n8:0",new double[]{6.5d},true,false));
    standards.get("d25:1").addAdduct("-H", 490.330299232d);
    standards.get("d25:1").addAdduct("H", 492.344851832d);
    standards.get("d25:1").addAdduct("-OH", 474.334287132d);
    standards.get("d25:1").addAdduct("Na", 514.326797532d);
    standards.put("m34:1", new ReferenceInfoVO("m18:1/n16:0",new double[]{20.4d},false,false));
    standards.get("m34:1").addAdduct("-H", 600.476235232d);
    standards.get("m34:1").addAdduct("H", 602.490787832d);
    standards.get("m34:1").addAdduct("-OH", 584.480223132d);
    standards.get("m34:1").addAdduct("Na", 624.472733532d);
    standards.put("d42:1", new ReferenceInfoVO("d18:1/n24:0",new double[]{30.4d},true,false));
    standards.get("d42:1").addAdduct("-H", 728.596350422d);
    standards.get("d42:1").addAdduct("H", 730.610903022d);
    standards.get("d42:1").addAdduct("-OH", 712.600338322d);
    standards.get("d42:1").addAdduct("Na", 752.592848722d);
    return standards;
  }

  
  public static LinkedHashMap<String,ReferenceInfoVO> getHexCerStandards(){
    LinkedHashMap<String,ReferenceInfoVO> standards = new LinkedHashMap<String,ReferenceInfoVO>();
    // the first boolean is if the standard shall be used in the evaluation, and the second if a position is detectable
    standards.put("d24:2", new ReferenceInfoVO("d18:2/n6:0",new double[]{5.4d},true,false));
    standards.get("d24:2").addAdduct("-H", 556.385491667d);
    standards.get("d24:2").addAdduct("HCOO", 602.390970997d);
    standards.get("d24:2").addAdduct("Cl", 592.362169967d);
    standards.get("d24:2").addAdduct("H", 558.400044267d);
    standards.get("d24:2").addAdduct("-OH", 540.389479567d);
    standards.get("d24:2").addAdduct("Na", 580.381989967d);
    standards.put("d36:2", new ReferenceInfoVO("d18:1/n18:1",new double[]{20.7d},true,false));
    standards.get("d36:2").addAdduct("-H", 724.573292507d);
    standards.get("d36:2").addAdduct("HCOO", 770.578771837d);
    standards.get("d36:2").addAdduct("Cl", 760.549970807d);
    standards.get("d36:2").addAdduct("H", 726.587845107d);
    standards.get("d36:2").addAdduct("-OH", 708.577280407d);
    standards.get("d36:2").addAdduct("Na", 748.569790807d);
    standards.put("d42:2", new ReferenceInfoVO("d18:1/n24:1",new double[]{26.6d},true,false));
    standards.get("d42:2").addAdduct("-H", 808.667192927d);
    standards.get("d42:2").addAdduct("HCOO", 854.672672257d);
    standards.get("d42:2").addAdduct("Cl", 844.643871227d);
    standards.get("d42:2").addAdduct("H", 810.681745527d);
    standards.get("d42:2").addAdduct("-OH", 792.671180827d);
    standards.get("d42:2").addAdduct("Na", 832.663691227d);
    standards.put("t44:0", new ReferenceInfoVO("t18:0/n26:0",new double[]{29.1d},true,false));
    standards.get("t44:0").addAdduct("-H", 856.724707837d);
    standards.get("t44:0").addAdduct("HCOO", 902.730187167d);
    standards.get("t44:0").addAdduct("Cl", 892.701386137d);
    standards.get("t44:0").addAdduct("H", 858.739260437d);
    standards.get("t44:0").addAdduct("-OH", 840.728695737d);
    standards.get("t44:0").addAdduct("Na", 880.721206137d);
    return standards;
  }
  
  public static LinkedHashMap<String,ReferenceInfoVO> getSMStandards(){
    LinkedHashMap<String,ReferenceInfoVO> standards = new LinkedHashMap<String,ReferenceInfoVO>();
    // the first boolean is if the standard shall be used in the evaluation, and the second if a position is detectable
    standards.put("d30:0", new ReferenceInfoVO("d18:0/n12:0",new double[]{14.5d},true,false));
    standards.get("d30:0").addAdduct("HCOO", 693.518828369d);
    standards.get("d30:0").addAdduct("Cl", 683.490027339d);
    standards.get("d30:0").addAdduct("H", 649.527901639d);
    standards.get("d30:0").addAdduct("Na", 671.509847339d);
    standards.put("d36:1", new ReferenceInfoVO("d18:1/n18:0",new double[]{21.2d},true,false));
    standards.get("d36:1").addAdduct("HCOO", 775.597078719d);
    standards.get("d36:1").addAdduct("Cl", 765.568277689d);
    standards.get("d36:1").addAdduct("H", 731.606151989d);
    standards.get("d36:1").addAdduct("Na", 753.588097689d);
    standards.put("d42:1", new ReferenceInfoVO("d18:1/n24:0",new double[]{27.6d},true,false));
    standards.get("d42:1").addAdduct("HCOO", 859.690979139d);
    standards.get("d42:1").addAdduct("Cl", 849.662178109d);
    standards.get("d42:1").addAdduct("H", 815.700052409d);
    standards.get("d42:1").addAdduct("Na", 837.681998109d);    
    return standards;
  }
 
}
