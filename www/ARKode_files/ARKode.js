// Created by iWeb 3.0.4 local-build-20120531

setTransparentGifURL('Media/transparent.gif');function applyEffects()
{var registry=IWCreateEffectRegistry();registry.registerEffects({shadow_2:new IWShadow({blurRadius:4,offset:new IWPoint(1.4142,1.4142),color:'#000000',opacity:0.500000}),shadow_1:new IWShadow({blurRadius:4,offset:new IWPoint(1.4142,1.4142),color:'#000000',opacity:0.500000}),shadow_0:new IWShadow({blurRadius:4,offset:new IWPoint(1.4142,1.4142),color:'#000000',opacity:0.500000}),stroke_0:new IWStrokeParts([{rect:new IWRect(-1,1,2,89),url:'ARKode_files/stroke.png'},{rect:new IWRect(-1,-1,2,2),url:'ARKode_files/stroke_1.png'},{rect:new IWRect(1,-1,257,2),url:'ARKode_files/stroke_2.png'},{rect:new IWRect(258,-1,2,2),url:'ARKode_files/stroke_3.png'},{rect:new IWRect(258,1,2,89),url:'ARKode_files/stroke_4.png'},{rect:new IWRect(258,90,2,3),url:'ARKode_files/stroke_5.png'},{rect:new IWRect(1,90,257,3),url:'ARKode_files/stroke_6.png'},{rect:new IWRect(-1,90,2,3),url:'ARKode_files/stroke_7.png'}],new IWSize(259,91))});registry.applyEffects();}
function hostedOnDM()
{return false;}
function onPageLoad()
{loadMozillaCSS('ARKode_files/ARKodeMoz.css')
adjustLineHeightIfTooBig('id1');adjustFontSizeIfTooBig('id1');adjustLineHeightIfTooBig('id2');adjustFontSizeIfTooBig('id2');adjustLineHeightIfTooBig('id3');adjustFontSizeIfTooBig('id3');adjustLineHeightIfTooBig('id4');adjustFontSizeIfTooBig('id4');Widget.onload();fixupAllIEPNGBGs();fixAllIEPNGs('Media/transparent.gif');applyEffects()}
function onPageUnload()
{Widget.onunload();}
