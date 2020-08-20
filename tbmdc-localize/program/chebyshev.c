#include  <math.h>
#include  <stdio.h>
#include  <stdlib.h>
#include  "md.h"
#include  "prototypes.h"

/* Chebyshev polynomials expansion coefficients up to 101 order */

void   get_chebyshev_coeff(void){
  double x, x2, x3, x4, x5;

  x  = ch.myu;  /* ch.myu = chemical potential */
  x2 = x*x;
  x3 = x2*x;
  x4 = x3*x;
  x5 = x4*x;

  /* [Hamiltonian] ----------------------------------------------------*/

      ch.H[1]=(0.997152953425486+0.099362199516628*x-2.628436812334222e-7*x2+0.0004053350587078263*x3-1.805444863907654e-6*x4+5.309838135171497e-6*x5);
      ch.H[2]=(-0.6336090037364808-7.5145124355548e-8*x+0.007755586336877767*x2+2.887400989895634e-6*x3+0.00004287025203956623*x4+2.782669686573931e-6*x5);
      ch.H[3]=(-0.002847006307478434-0.0992600402288005*x-3.645807424007939e-7*x2+0.001212668052181541*x3-2.506490853926179e-6*x4+8.60259004813889e-6*x5);
      ch.H[4]=(0.2145631535509881-1.679654523066751e-7*x-0.02322165251896424*x2+6.43397959384319e-6*x3+0.0002262160549422853*x4+6.326093163901695e-6*x5);
      ch.H[5]=(-0.002846876705346845+0.0989540364169679*x-1.056677072705578e-6*x2-0.006030441213931537*x3-7.311571861399357e-6*x4+0.00005516639496393176*x5);
      ch.H[6]=(-0.1236618798585157-1.390474873622747e-6*x+0.03853227136642338*x2+0.00005297344750695696*x3-0.001731519946997195*x4+0.00005407185255809111*x5);
      ch.H[7]=(-0.002846461533252191-0.0984487639672454*x-0.00001097346124731672*x2+0.01406548188174787*x3-0.0000767567721077622*x4-0.000425223931968753*x5);
      ch.H[8]=(0.0926540407195333+0.00001341306248337098*x-0.05353838500978458*x2-0.0005060557651791325*x3+0.005736572095128973*x4-0.0005499455172636081*x5);
      ch.H[9]=(-0.002849435409407233+0.0977606965903382*x+0.0001314662704700413*x2-0.02569822443282403*x3+0.000934314178621447*x4+0.001225917694173622*x5);
      ch.H[10]=(-0.0664297417713815-0.00004814375540652608*x+0.06808810572941172*x2+0.001788134660338529*x3-0.01375320287157576*x4+0.00214258150283126*x5);

      ch.H[11]=(-0.002832652167940699-0.0969262736386177*x-0.0005906849088626514*x2+0.04198274273906294*x3-0.004290128454561124*x4-0.002207237325522289*x5);
      ch.H[12]=(0.05894300727897226+0.0001109986219474054*x-0.082123270895594*x2-0.004015770039566216*x3+0.02699659173023896*x4-0.005608992110704916*x5);
      ch.H[13]=(-0.00288576061893502+0.0960048196698706*x+0.001756782083586701*x2-0.06473696540648394*x3+0.01313153079974643*x4+0.002555617301397824*x5);
      ch.H[14]=(-0.04403973089786133-0.0001849934094800096*x+0.0958096065362439*x2+0.006365725893498172*x3-0.04554666208381347*x4+0.0115054297669659*x5);
      ch.H[15]=(-0.00275257605788965-0.0950719810051539*x-0.004050549384107265*x2+0.0962411184103548*x3-0.03143490991430999*x4-0.000892716618039209*x5);
      ch.H[16]=(0.04291589069916691+0.0002153610217232091*x-0.1096468603084503*x2-0.006490887984820437*x3+0.06749169807215382*x4-0.01977693790209641*x5);
      ch.H[17]=(-0.003022851035915565+0.0942023145951485*x+0.007768750218236465*x2-0.1385028456165474*x3+0.06331948523008358*x4-0.004486878243665567*x5);
      ch.H[18]=(-0.03199361290131086-0.0001010694630733407*x+0.1244965664128236*x2+0.0002515461700852181*x3-0.0884039445882429*x4+0.02945188789284553*x5);
      ch.H[19]=(-0.002547718650653497-0.0934436128573671*x-0.0128515049132701*x2+0.1922153879205812*x3-0.1117511930668034*x4+0.01516688250176742*x5);
      ch.H[20]=(0.03357151838181691-0.0003009551228886428*x-0.1414797745738603*x2+0.01800813412333101*x3+0.1014536541406827*x4-0.03855232762527719*x5);

      ch.H[21]=(-0.00327661537104958+0.0927891946678324*x+0.01866225761279582*x2-0.2556849480823788*x3+0.1771096059900122*x4-0.03205422795354868*x5);
      ch.H[22]=(-0.02462483762032752+0.001156169482521655*x+0.1617327042255261*x2-0.0545485526942164*x3-0.0983123488395739*x4+0.04429893806549847*x5);
      ch.H[23]=(-0.002281802344257454-0.0921564871848134*x-0.02386330804097867*x2+0.324081166244115*x3-0.256041590620427*x4+0.05487274271692427*x5);
      ch.H[24]=(0.0277230465145002-0.002619022247655099*x-0.1860550197426895*x2+0.114791634535978*x3+0.0707694156889818*x4-0.04360741426169172*x5);
      ch.H[25]=(-0.003476321508366489+0.0913799449716343*x+0.02645839360432088*x2-0.3893370306044479*x3+0.3410144214576468*x4-0.0818672960494679*x5);
      ch.H[26]=(-0.01999688987118397+0.004792054488251228*x+0.2145288852911136*x2-0.2016393968932248*x3-0.0127544871272916*x4+0.03378468862309346*x5);
      ch.H[27]=(-0.002237652321927016-0.090223576462657*x-0.02403698761351358*x2+0.4408897486456774*x3-0.420841835193514*x4+0.1098386851988361*x5);
      ch.H[28]=(0.02406565275207701-0.007686947194227501*x-0.2462127121154347*x2+0.3139727822571374*x3-0.07770806077228331*x4-0.01326467278170536*x5);
      ch.H[29]=(-0.003252789981193588+0.0884136990478985*x+0.01420161369119738*x2-0.4672435798117147*x3+0.4822213162126046*x4-0.1345499027575248*x5);
      ch.H[30]=(-0.01709139109556501+0.01119786441487027*x+0.279011618521776*x2-0.4457673391837321*x3+0.1971328307657282*x4-0.01780498639411755*x5);

      ch.H[31]=(-0.002813286377756013-0.0856872280946706*x+0.004895919168095653*x2+0.4581054668052989*x3-0.5120517862314503*x4+0.1514477156406239*x5);
      ch.H[32]=(0.02167230221438468-0.01509557922007231*x-0.3097948882263775*x2+0.5861586472242978*x3-0.3362219022365874*x4+0.05729161423218753*x5);
      ch.H[33]=(-0.002260375777745774+0.081846265711702*x-0.03402480678898461*x2-0.406659137927467*x3+0.5000718656135667*x4-0.1565560759654415*x5);
      ch.H[34]=(-0.01497967460371848+0.01904660151215728*x+0.3347756286345859*x2-0.7205857304712905*x3+0.4807388152522769*x4-0.101220264379515*x5);
      ch.H[35]=(-0.004213145616994423-0.07680726493204761*x+0.07248191113117565*x2+0.3114576946201483*x3-0.4412348486750362*x4+0.1473438599028096*x5);
      ch.H[36]=(0.01959871996162825-0.02265587310252021*x-0.3501029386650715*x2+0.832885001302867*x3-0.613533679528315*x4+0.1443323643612712*x5);
      ch.H[37]=(-0.0005117167072591747+0.07063335083048226*x-0.1179461654692766*x2-0.1774572749709132*x3+0.3372549336492414*x4-0.1233639388471902*x5);
      ch.H[38]=(-0.01274819472103392+0.02552597109800369*x+0.3525590943342135*x2-0.90796603482883*x3+0.717363383866783*x4-0.1809613510899701*x5);
      ch.H[39]=(-0.006164674219241891-0.06354149836191745*x+0.166643641096342*x2+0.01587988371967384*x3-0.196919792368716*x4+0.0865058164807444*x5);
      ch.H[40]=(0.01704685830893888-0.02732159049420447*x-0.3402191484742821*x2+0.934536861110495*x3-0.7779451838706213*x4+0.206047506462726*x5);

      ch.H[41]=(0.001465217332094962+0.05588148390860794*x-0.2138187945073798*x2+0.1571561087087635*x3+0.03502959817600758*x4-0.04078890766203347*x5);
      ch.H[42]=(-0.00981162639733959+0.02782641919196069*x+0.3129279460403494*x2-0.907304263784143*x3+0.786605290831881*x4-0.2160771068795354*x5);
      ch.H[43]=(-0.00794083646350616-0.04808959631742673*x+0.2544373495375396*x2-0.3230244452847372*x3+0.1298724359258757*x4-0.00827316620771699*x5);
      ch.H[44]=(0.01372678323076661-0.02698079606760825*x-0.2724816269177758*x2+0.828166743556049*x3-0.7419564472342006*x4+0.2097446858903233*x5);
      ch.H[45]=(0.002845428160451121+0.04062556338982364*x-0.2839907628520262*x2+0.4635526947145719*x3-0.2785423376824551*x4+0.05461056037257248*x5);
      ch.H[46]=(-0.006220264487552028+0.02489244878734811*x+0.2224602480103368*x2-0.7061234883113129*x3+0.6502422959621375*x4-0.1881977472956848*x5);
      ch.H[47]=(-0.00874694544326694-0.03390473981260177*x+0.2992449764699182*x2-0.5639775478277606*x3+0.3942431263821601*x4-0.0926444580836414*x5);
      ch.H[48]=(0.01002126053094116-0.02181820502233773*x-0.1677313678556624*x2+0.5558931756640562*x3-0.5242720336750266*x4+0.1548152895639322*x5);
      ch.H[49]=(0.002995528446019309+0.02823851252375376*x-0.2987839152864284*x2+0.6153098843659916*x3-0.4654770299382348*x4+0.1182292858475363*x5);
      ch.H[50]=(-0.002625871143056866+0.01812054297736353*x+0.1137141899496822*x2-0.3955080487138808*x3+0.3811682055824418*x4-0.1145725297070798*x5);

      ch.H[51]=(-0.00821868906260691-0.02379396847478387*x+0.2832391560742881*x2-0.6156507843589724*x3+0.4877103125427527*x4-0.1293028735017621*x5);
      ch.H[52]=(0.006741075120240487-0.01420783640189778*x-0.0655440832426626*x2+0.2433589901796012*x3-0.2393942930135513*x4+0.0731305260442237*x5);
      ch.H[53]=(0.00188317105592187+0.02057966006914855*x-0.2551617188440526*x2+0.5702226117381242*x3-0.4637815830488721*x4+0.1261252575607284*x5);
      ch.H[54]=(0.0001387848059766208+0.01047008190096019*x+0.02729684678564203*x2-0.1152701485289336*x3+0.1156634828118789*x4-0.03584186797134398*x5);
      ch.H[55]=(-0.006661085700920922-0.01845882303646119*x+0.2185650620958988*x2-0.4901471931019921*x3+0.402966198535427*x4-0.1110777347066642*x5);
      ch.H[56]=(0.004602222474366676-0.007222222043734106*x-0.001417255958309226*x2+0.02215759583989212*x3-0.02232998599005443*x4+0.006872019990137979*x5);
      ch.H[57]=(0.0000977591302591643+0.01718591308947737*x-0.1782336243239251*x2+0.3902547270580754*x3-0.3189445956079183*x4+0.0880862052774889*x5);
      ch.H[58]=(0.001589304620801982+0.004664966782792784*x-0.01154846660596604*x2+0.03131395836428762*x3-0.03426390806787138*x4+0.01139914926837201*x5);
      ch.H[59]=(-0.004847436896001319-0.01645802120310803*x+0.1389345892908216*x2-0.2863877915502924*x3+0.2271306201199258*x4-0.06180612215463013*x5);
      ch.H[60]=(0.003787130090498194-0.002868902169222231*x+0.01288062081881171*x2-0.04688282742770681*x3+0.05423282825642927*x4-0.01859916758447943*x5);

      ch.H[61]=(-0.001531474255332112+0.01597045060529812*x-0.1046834791208261*x2+0.1927390164791152*x3-0.1419155793708196*x4+0.03674555090477742*x5);
      ch.H[62]=(0.001859995491840417+0.001782672024434556*x-0.005307218874708398*x2+0.03170499960993914*x3-0.04339290170990187*x4+0.01619540098328824*x5);
      ch.H[63]=(-0.003532200620785042-0.01546582933012414*x+0.07819543326404311*x2-0.1197189656770942*x3+0.07436034957195731*x4-0.01650180858056562*x5);
      ch.H[64]=(0.003885484268510136-0.00126125761462831*x-0.007590889162664304*x2+0.003492806642596526*x3+0.01155304965224349*x4-0.007005512786832241*x5);
      ch.H[65]=(-0.002434182600054432+0.01476838194763748*x-0.06060942090221765*x2+0.0727092365829371*x3-0.030735417968611*x4+0.003249235899566213*x5);
      ch.H[66]=(0.001551860369278233+0.001106849844586752*x+0.02206444954024596*x2-0.0466981366765224*x3+0.02973043461185906*x4-0.005487438148400027*x5);
      ch.H[67]=(-0.003032274961948873-0.013798730772817*x+0.05151326305826261*x2-0.05185213420946159*x3+0.01210423076623537*x4+0.002447839259507039*x5);
      ch.H[68]=(0.004216673823440382-0.001113141804791714*x-0.03483570784502538*x2+0.0868627598743019*x3-0.06945071777013422*x4+0.01786121437931934*x5);
      ch.H[69]=(-0.002566607563753515+0.01256891282159442*x-0.0492367465392416*x2+0.0528105766499917*x3-0.01491842330094762*x4-0.001496004163511403*x5);
      ch.H[70]=(0.001310757889089102+0.001104238117268083*x+0.04357980566184727*x2-0.1156959526913689*x3+0.0990511544572566*x4-0.02737526956780504*x5);

      ch.H[71]=(-0.003147607148180979-0.01116122152692049*x+0.05133107335742427*x2-0.06824727494884007*x3+0.03239468439205033*x4-0.004127290793289978*x5);
      ch.H[72]=(0.00426593055466607-0.000961496872687705*x-0.04717508264722763*x2+0.1287121557125205*x3-0.1136087986901586*x4+0.03237684732916964*x5);
      ch.H[73]=(-0.002304310332234594+0.00969723037661121*x-0.05512466540141939*x2+0.0896565267765293*x3-0.05631165350733606*x4+0.01192389627718418*x5);
      ch.H[74]=(0.001458040811309487+0.0006348579272535593*x+0.04570190138518751*x2-0.1254225340798279*x3+0.1121591235537326*x4-0.03243898000532919*x5);
      ch.H[75]=(-0.003417786063255428-0.00830448237603913*x+0.05824194143582284*x2-0.109152820328617*x3+0.07882238353446462*x4-0.01944682884100096*x5);
      ch.H[76]=(0.00392574968304881-0.0001387634023856074*x-0.04022080252481224*x2+0.1087422936022532*x3-0.097209599850523*x4+0.02823802548822751*x5);
      ch.H[77]=(-0.002098863400371091+0.007087784628393844*x-0.05899145345988422*x2+0.1208752021287752*x3-0.0939232642063464*x4+0.0247766530475588*x5);
      ch.H[78]=(0.001907826355223182-0.000464139600652269*x+0.0324004637743849*x2-0.0838331961960779*x3+0.07364767256294069*x4-0.02123067451079213*x5);
      ch.H[79]=(-0.003478948302095434-0.006110145513550964*x+0.05656629809500074*x2-0.1217846570614897*x3+0.0983346187634584*x4-0.02681958512213493*x5);
      ch.H[80]=(0.003422901001260842+0.001086298881855072*x-0.02408704869661659*x2+0.05668645007661753*x3-0.04734538007066693*x4+0.01322382313450216*x5);

      ch.H[81]=(-0.002175370753226261+0.005385700961185978*x-0.05104368535460817*x2+0.1117817953039642*x3-0.0916994471492895*x4+0.02539203579178794*x5);
      ch.H[82]=(0.002362949066299996-0.001639671489974281*x+0.01690671232744256*x2-0.03276404778450219*x3+0.02378504833732971*x4-0.005940427613182532*x5);
      ch.H[83]=(-0.003261967914281837-0.004884164334199705*x+0.04321203503458733*x2-0.0932188725166184*x3+0.07615866271977182*x4-0.02110523250630935*x5);
      ch.H[84]=(0.003052372411081037+0.002056010685502254*x-0.01197458361308168*x2+0.01596249131583317*x3-0.006982393745001571*x4+0.0006712612959389335*x5);
      ch.H[85]=(-0.00246225461094778+0.004544012089495257*x-0.0342846673129655*x2+0.06999321841746883*x3-0.05548144199666767*x4+0.01510326889922195*x5);
      ch.H[86]=(0.002590962954138157-0.002299436986453708*x+0.00975245421217157*x2-0.00806111490559602*x3-0.00111616642922333*x4+0.001927932929594258*x5);
      ch.H[87]=(-0.002938973663073199-0.0042902116348821*x+0.0255739645671932*x2-0.04646889895678314*x3+0.03399643551095212*x4-0.00873109742856035*x5);
      ch.H[88]=(0.002945954253513049+0.002369876814017289*x-0.0100616394425364*x2+0.0086951735829071*x3+0.0007085473812135941*x4-0.001873332925647143*x5);
      ch.H[89]=(-0.002740345673762576+0.004051982342037232*x-0.01819766020574135*x2+0.02647401285419894*x3-0.0155761176030664*x4+0.003212789654355247*x5);
      ch.H[90]=(0.002572374401671036-0.002298107646564531*x+0.01222505699019486*x2-0.01577979675988733*x3+0.006251212173525861*x4-0.0002586185380061126*x5);

      ch.H[91]=(-0.002721964858260818-0.003776770593377654*x+0.01287224989347919*x2-0.01256991622009261*x3+0.002879019587493567*x4+0.0005934191229014401*x5);
      ch.H[92]=(0.003024027347763662+0.002134827848376548*x-0.01528983002988517*x2+0.02622785751159208*x3-0.01675200299410836*x4+0.003545025030721191*x5);
      ch.H[93]=(-0.002848725533511187+0.003438024291889839*x-0.00982274985294395*x2+0.005703474490583507*x3+0.00302965807727488*x4-0.002322864652177325*x5);
      ch.H[94]=(0.002454381012731703-0.001937079184836821*x+0.01827216421317033*x2-0.0367662462688709*x3+0.02753981629316195*x4-0.006971149428486689*x5);
      ch.H[95]=(-0.002694542105210992-0.003036065927838565*x+0.0088097352572848*x2-0.005256783674788028*x3-0.002656522532921681*x4+0.00210557983809094*x5);
      ch.H[96]=(0.003112445705866995+0.001755405165172216*x-0.02037013409935506*x2+0.04466314332803814*x3-0.03585581567216192*x4+0.0096648298067323*x5);
      ch.H[97]=(-0.002787019173071652+0.002592970208365638*x-0.00925005781529599*x2+0.00942279440211575*x3-0.002276201628211036*x4-0.0004600353540348083*x5);
      ch.H[98]=(0.002398162939701786-0.001624441559252099*x+0.02110471727078769*x2-0.04822553426139574*x3+0.0399518431745848*x4-0.01106315786062166*x5);
      ch.H[99]=(-0.002786384871161688-0.002143518880742416*x+0.01039152485106502*x2-0.01577678459128056*x3+0.00939351321228591*x4-0.001880970193332432*x5);
      ch.H[100]=(0.003094008765667559+0.001558465471887078*x-0.02037035225732303*x2+0.04699608694107496*x3-0.03930660026602348*x4+0.01098566961542682*x5);
      ch.H[101]=(-0.002670751032911042+0.001724842204585286*x-0.01149637037420964*x2+0.0218913240524708*x3-0.01628386144084718*x4+0.004163962598525191*x5);

  /* [Force] ----------------------------------------------------------*/


      ch.F[1]=(0.997152953166889+0.0993109700764858*x-2.49969207922139e-7*x2+0.0004034077162742047*x3-1.749312145876739e-6*x4+5.227019758690783e-6*x5);
      ch.F[2]=(-0.6339366079472219-7.378761006389184e-8*x+0.00774358616609247*x2+2.820991005086637e-6*x3+0.00004259485540324933*x4+2.728101042955994e-6*x5);
      ch.F[3]=(-0.002847006499970031-0.0994135057823221*x-3.560837925965834e-7*x2+0.001209493320590613*x3-2.446131875245523e-6*x4+8.49273761668224e-6*x5);
      ch.F[4]=(0.2155431314811472-1.652229764728746e-7*x-0.02328154365984975*x2+6.330531985029713e-6*x3+0.000225518788752502*x4+6.229291610907569e-6*x5);
      ch.F[5]=(-0.002846877152985215+0.0997180533757178*x-1.03695192775977e-6*x2-0.006052366325151846*x3-7.17703964457107e-6*x4+0.00005486104402579212*x5);
      ch.F[6]=(-0.1252857905904949-1.381356266949506e-6*x+0.0389492868880813*x2+0.0000526219870575413*x3-0.001738644072238228*x4+0.00005373553223324845*x5);
      ch.F[7]=(-0.002846463372286409-0.1002186493467138*x-0.00001089246309987858*x2+0.01426012126305558*x3-0.0000762004723816474*x4-0.0004281251156015843*x5);
      ch.F[8]=(0.0949079737646549+0.00001351779128738892*x-0.05477773222925398*x2-0.0005100512115645392*x3+0.005826530068817088*x4-0.0005539813639828038*x5);
      ch.F[9]=(-0.002849460754896981+0.1009147559412115*x+0.000132581691363971*x2-0.02641171483965239*x3+0.000942098078925332*x4+0.001256054440531572*x5);
      ch.F[10]=(-0.06929446120849004-0.00004945371352051835*x+0.07077647353472309*x2+0.001837459845097367*x3-0.01419227297311351*x4+0.002196844195201802*x5);

      ch.F[11]=(-0.002832265812821061-0.1018211582363161*x-0.0006076689679754201*x2+0.04388503615555839*x3-0.004411140119144917*x4-0.002315221752034472*x5);
      ch.F[12]=(0.0623941863341857+0.0001167812026542323*x-0.0870260969738447*x2-0.004228873730553704*x3+0.02842511571045613*x4-0.005876345342547542*x5);
      ch.F[13]=(-0.002887896040346758+0.1029744341798566*x+0.001850508976143468*x2-0.06904741865738045*x3+0.01381796764909129*x4+0.002773177640115805*x5);
      ch.F[14]=(-0.04804864267002879-0.0002005135470567339*x+0.1038218591545225*x2+0.006916113139251304*x3-0.04912701456080162*x4+0.01236002558977209*x5);
      ch.F[15]=(-0.002744973953096686-0.1044309301490884*x-0.004383497525723193*x2+0.1050783565242176*x3-0.0339624119304317*x4-0.0010962950268877*x5);
      ch.F[16]=(0.04745115139744556+0.0002430618105443179*x-0.1218254865115581*x2-0.007391018322577882*x3+0.07484583269564768*x4-0.02185536458629353*x5);
      ch.F[17]=(-0.003043352039026881+0.1062527175789484*x+0.00866390379117056*x2-0.155293929385822*x3+0.07043997726155315*x4-0.004808267946437522*x5);
      ch.F[18]=(-0.03702451726157601-0.0001275426055279648*x+0.1421601819106024*x2+0.000808891826430024*x3-0.1011394320323047*x4+0.03357747173008369*x5);
      ch.F[19]=(-0.002502664552731493-0.1084811020745276*x-0.01481013226539388*x2+0.2219803961582969*x3-0.128319723251417*x4+0.01718541660493144*x5);
      ch.F[20]=(0.03907300934292052-0.0003280549759063927*x-0.1663787676457536*x2+0.02011530475512847*x3+0.1201561075045132*x4-0.04546463163604677*x5);

      ch.F[21]=(-0.00336036628521148+0.1111013566242358*x+0.0222792612541337*x2-0.3049208201379981*x3+0.2103464898070356*x4-0.03780472516786454*x5);
      ch.F[22]=(-0.03058325898133289+0.001359273331134523*x+0.196252746094255*x2-0.06446222666459927*x3-0.1210713527008745*x4+0.0541806931055548*x5);
      ch.F[23]=(-0.002147910719386286-0.1140064433420635*x-0.02958546890764869*x2+0.3999670449445289*x3-0.3150533901937957*x4+0.06723781175353859*x5);
      ch.F[24]=(0.03414082109941198-0.003222891368089984*x-0.2333716417923885*x2+0.1416152135025923*x3+0.0915561120597026*x4-0.0554918116589079*x5);
      ch.F[25]=(-0.003660075866840958+0.1169702919347577*x+0.03416867497238058*x2-0.4981132450861532*x3+0.4353141769488199*x4-0.1042262644130854*x5);
      ch.F[26]=(-0.02689359111227353+0.006145868872404682*x+0.2786002065627759*x2-0.2590521773208076*x3-0.02018975129075741*x4+0.04501551947510289*x5);
      ch.F[27]=(-0.002027826457649028-0.1196415250774628*x-0.0325127821219525*x2+0.5855164731024694*x3-0.5578828075411548*x4+0.1453446087878396*x5);
      ch.F[28]=(0.03147303499520559-0.01026343379332188*x-0.3315013239061905*x2+0.4198097991630539*x3-0.100569145690333*x4-0.01923168824037374*x5);
      ch.F[29]=(-0.003429524079531267+0.1215656845763876*x+0.02056865079828346*x2-0.6447594133196417*x3+0.6643252858239718*x4-0.1851251171507394*x5);
      ch.F[30]=(-0.02504105173334285+0.01555929084811227*x+0.3898758462240877*x2-0.6202090265399078*x3+0.2714953652852084*x4-0.02345326592801672*x5);

      ch.F[31]=(-0.002772162889469601-0.1222381944163051*x+0.005536032280234328*x2+0.6573758226215917*x3-0.7334564450542587*x4+0.2167177629018233*x5);
      ch.F[32]=(0.03017665795460127-0.02182341316116503*x-0.4495794567638407*x2+0.848505848095775*x3-0.4845565663967425*x4+0.0818543780935148*x5);
      ch.F[33]=(-0.002020097469694528+0.1211826352241047*x-0.0486034101619185*x2-0.607317430307453*x3+0.7450105744529428*x4-0.2330273590087679*x5);
      ch.F[34]=(-0.02400936082939949+0.02864155922772806*x+0.5047406939425879*x2-1.085031952589764*x3+0.7224130065164611*x4-0.1516255383951923*x5);
      ch.F[35]=(-0.004912649877948357-0.1180412932120673*x+0.1094669389080603*x2+0.4847212862500583*x3-0.6838965808307198*x4+0.228133270510898*x5);
      ch.F[36]=(0.02906162243427952-0.03542527589444932*x-0.5484253523029451*x2+1.304113870734732*x3-0.959876651258113*x4+0.2255229229273628*x5);
      ch.F[37]=(0.000833788126257426+0.1126595512061555*x-0.1863060496967917*x2-0.2891354145623342*x3+0.544125536324865*x4-0.1986935150902495*x5);
      ch.F[38]=(-0.02247632098726821+0.04148260797788239*x+0.573685891687408*x2-1.477657255342501*x3+1.167318792505617*x4-0.2943514195636335*x5);
      ch.F[39]=(-0.00831948857844464-0.1051443158733429*x+0.2743801025815904*x2+0.03135344522538614*x3-0.3314411904904398*x4+0.1449978784638712*x5);
      ch.F[40]=(0.02679661063179743-0.04611996558621077*x-0.5748261248978273*x2+1.579849881727973*x3-1.315529133672353*x4+0.3484643673209895*x5);

      ch.F[41]=(0.004533083587412059+0.0958801871180497*x-0.3663187375215011*x2+0.2667761057602596*x3+0.06387328558209032*x4-0.07138142731222965*x5);
      ch.F[42]=(-0.01928065456021291+0.04875708465398437*x+0.5486311576585894*x2-1.59208785163467*x3+1.381113666260565*x4-0.3795361999066223*x5);
      ch.F[43]=(-0.01193391846023604-0.0854950522178224*x+0.4529912957495294*x2-0.5747952128511371*x3+0.2301557913054222*x4-0.01416425515448051*x5);
      ch.F[44]=(0.02258690898182824-0.04903211611584525*x-0.4952856399653413*x2+1.507067601807845*x3-1.35129209331654*x4+0.3822347427980061*x5);
      ch.F[45]=(0.007664504166199348+0.07477754631137539*x-0.5248435212539539*x2+0.85805907954439*x3-0.5159974095714996*x4+0.1011668933068928*x5);
      ch.F[46]=(-0.01416185967935656+0.04687415175589741*x+0.4187428334334653*x2-1.331072125747091*x3+1.226993606090155*x4-0.3554137814458051*x5);
      ch.F[47]=(-0.01418015751281196-0.06455981436305203*x+0.5734685647371419*x2-1.083374291261279*x3+0.7584987717855137*x4-0.1784345314857205*x5);
      ch.F[48]=(0.01680102276470633-0.04252621230779173*x-0.3264087698483747*x2+1.083810951039557*x3-1.023461118598964*x4+0.3025339381008781*x5);
      ch.F[49]=(0.0087384225517088+0.05558748709969973*x-0.5931037699553431*x2+1.224753709392429*x3-0.928158414912598*x4+0.2360521576186841*x5);     
      ch.F[50]=(-0.00810817357904828+0.03651174275998293*x+0.2281488431373706*x2-0.7956818978995309*x3+0.76809480064288*x4-0.2311829725662332*x5);

      ch.F[51]=(-0.01391360584506819-0.04840253238027488*x+0.5817372274411402*x2-1.268094413212738*x3+1.006436218983181*x4-0.26719479915922*x5);
      ch.F[52]=(0.01092358033551921-0.02954968635758028*x-0.1347753880912177*x2+0.5028906611420902*x3-0.4958745696925447*x4+0.1517669441452929*x5);
      ch.F[53]=(0.007171454628125574+0.04326259566686346*x-0.5415745095574822*x2+1.213766631733348*x3-0.989079051209418*x4+0.2693663729826431*x5);
      ch.F[54]=(-0.002878263797493542+0.02243438154625091*x+0.05629730704344408*x2-0.2413542340331519*x3+0.243250425587836*x4-0.07563230183181101*x5);
      ch.F[55]=(-0.01123824717849647-0.04011297295211767*x+0.4787435342152278*x2-1.076525129728463*x3+0.886733501998144*x4-0.2447948979716778*x5);
      ch.F[56]=(0.006704050803359948-0.01590415943730832*x-0.0002757097943562104*x2+0.0406047457338222*x3-0.0417438309728007*x4+0.01300935948841936*x5);
      ch.F[57]=(0.003760560381593625+0.03861619396286639*x-0.4022755106657179*x2+0.88272338816053*x3-0.7227026173467645*x4+0.1999001693925749*x5);
      ch.F[58]=(0.00007825261674579209+0.01052479467653075*x-0.02938394564231023*x2+0.0810719694923621*x3-0.0874248199336826*x4+0.02879220821904917*x5);
      ch.F[59]=(-0.007523216645842351-0.03823196582659377*x+0.3225551001893502*x2-0.6653913436721831*x3+0.5283247765663431*x4-0.1439548278057343*x5);
      ch.F[60]=(0.005051352342391401-0.006610383001589368*x+0.03296077005483833*x2-0.1192902132137958*x3+0.1365957049702151*x4-0.04656450573585636*x5);

      ch.F[61]=(0.0002207083636021227+0.03832979995355741*x-0.2495467525679041*x2+0.4581925166061308*x3-0.3369623649073453*x4+0.0872347533824255*x5);
      ch.F[62]=(0.0005360759416447208+0.004195711511634711*x-0.01466859289617374*x2+0.0835166354386907*x3-0.1122008653617642*x4+0.04153117659201242*x5);
      ch.F[63]=(-0.004540972583327036-0.03831038999501643*x+0.1911500835308588*x2-0.2895005463322528*x3+0.1778672708292269*x4-0.03905891126655411*x5);
      ch.F[64]=(0.005490359246712141-0.003062838473440281*x-0.01830838150487546*x2+0.005590747546557865*x3+0.03255298821528174*x4-0.01878948315402935*x5);
      ch.F[65]=(-0.001911761098443029+0.03771103670813464*x-0.152005022510732*x2+0.1777800690454512*x3-0.07117919233837366*x4+0.006129470189626891*x5);
      ch.F[66]=(-0.0004404302180731709+0.002813102578214413*x+0.05728239688238224*x2-0.1208185235802037*x3+0.07653537433592008*x4-0.01397782539020223*x5);
      ch.F[67]=(-0.003350376389822513-0.03627509229576426*x+0.1329658244122928*x2-0.1291417599372356*x3+0.02502222365801735*x4+0.00849918695631345*x5);
      ch.F[68]=(0.006582009843236074-0.002966759644953262*x-0.0937213055258351*x2+0.2339785316562937*x3-0.1872906049632971*x4+0.04821272962093134*x5);
      ch.F[69]=(-0.002198890644032429+0.03397420921513698*x-0.1313165765041662*x2+0.1374463092054106*x3-0.03517218802201016*x4-0.005646236420093076*x5);
      ch.F[70]=(-0.001309943548528133+0.003067818738127357*x+0.1207621756125119*x2-0.3212520894482327*x3+0.2754679504377349*x4-0.07623153234902644*x5);

      ch.F[71]=(-0.003747250312680185-0.03098291283388416*x+0.1416439151482884*x2-0.186776417769568*x3+0.0871932160470942*x4-0.01056739653863123*x5);
      ch.F[72]=(0.006957172386641845-0.00277221444974138*x-0.1342878103234913*x2+0.3672311289694679*x3-0.3246667684586481*x4+0.0926459606565946*x5);
      ch.F[73]=(-0.001387857222785754+0.02761513806076046*x-0.1571508101718401*x2+0.2556078008044939*x3-0.1604218167675615*x4+0.03390910747670486*x5);
      ch.F[74]=(-0.001097988204933059+0.001902751947476635*x+0.1333815949739908*x2-0.3669303913154253*x3+0.3286696531196749*x4-0.0951837408503048*x5);
      ch.F[75]=(-0.004634435329975618-0.02423964089516675*x+0.1711166850250058*x2-0.3216958742737466*x3+0.2328250688886519*x4-0.05754273884033411*x5);
      ch.F[76]=(0.006165939101971254-0.0004627769154524848*x-0.1201168734083503*x2+0.3255230195552783*x3-0.2914838598597462*x4+0.0847865232416048*x5);
      ch.F[77]=(-0.0006833112403477599+0.02119418929070345*x-0.1781984457121185*x2+0.3666145347678804*x3-0.2856455584852895*x4+0.07550875874957354*x5);
      ch.F[78]=(0.0001025946466152362-0.001389899247199909*x+0.0987786051639665*x2-0.2560520133581642*x3+0.225280424039709*x4-0.06502925192950469*x5);
      ch.F[79]=(-0.004939915685834473-0.01871681344511453*x+0.1753225669348838*x2-0.3790401511924201*x3+0.3068789722817506*x4-0.0838649344721652*x5);
      ch.F[80]=(0.004755506877898899+0.003398882121387583*x-0.07472605995426184*x2+0.1757789678075609*x3-0.1468614586481891*x4+0.04104866225541034*x5);

      ch.F[81]=(-0.000834690821858109+0.01690689876152932*x-0.1620236955976199*x2+0.356216161766184*x3-0.2929629142267257*x4+0.0812768207309853*x5);
      ch.F[82]=(0.001458417496177444-0.005279500841962041*x+0.05316334591230887*x2-0.1020647272180955*x3+0.07354764088223224*x4-0.0182587237062511*x5);
      ch.F[83]=(-0.004342790009617075-0.01572110642944479*x+0.1402136063496544*x2-0.3034714610315059*x3+0.2484974087919979*x4-0.06898656601611024*x5);
      ch.F[84]=(0.003638809686959987+0.006786626458730172*x-0.03808069803230652*x2+0.04869561974426782*x3-0.01942091275193047*x4+0.001151461457668019*x5);
      ch.F[85]=(-0.001726020151676396+0.01500096809097353*x-0.1134858147672882*x2+0.2320941234116265*x3-0.1842606068473443*x4+0.0502313290510211*x5);
      ch.F[86]=(0.002176301909198671-0.007765460197074868*x+0.03156929647832595*x2-0.02338694032692788*x3-0.007309470406209644*x4+0.007566051191348128*x5);
      ch.F[87]=(-0.003314582413649156-0.01452238589043728*x+0.0861495500917901*x2-0.1561864279204375*x3+0.1141278939910823*x4-0.02929636307912027*x5);
      ch.F[88]=(0.003337561722471687+0.00817610845317932*x-0.03361375934890284*x2+0.02686993668734427*x3+0.005370407869680358*x4-0.007346779449143284*x5);
      ch.F[89]=(-0.002675229954882812+0.01405351795680911*x-0.06222544035111168*x2+0.0892845251504258*x3-0.05168036521485315*x4+0.01044946458252421*x5);
      ch.F[90]=(0.002083649647904049-0.00808960273095343*x+0.04235434063122522*x2-0.05358486704182661*x3+0.02014012660525845*x4-0.0003375660152122375*x5);

      ch.F[91]=(-0.002576196313485509-0.01340719887943504*x+0.04462002050701786*x2-0.04149656533362072*x3+0.007275046702339583*x4+0.002996050352909241*x5);
      ch.F[92]=(0.003679934108806404+0.007659649397103365*x-0.0547111363311174*x2+0.0936515288977827*x3-0.05964139190271257*x4+0.01256533888357708*x5);
      ch.F[93]=(-0.003082072428639667+0.01247685519903437*x-0.03463529721139081*x2+0.01772080654614957*x3+0.0138343425448716*x4-0.00929756880302476*x5);
      ch.F[94]=(0.001619749349302634-0.007079329454271324*x+0.06719606710069017*x2-0.1355451037160858*x3+0.1017175055564789*x4-0.02578703952743274*x5);
      ch.F[95]=(-0.002492339750005621-0.01124985665073353*x+0.03188081718185657*x2-0.01720881734159052*x3-0.01203284902755569*x4+0.00847861831063793*x5);
      ch.F[96]=(0.004084020609877018+0.006534274607023955*x-0.07671655028687681*x2+0.16881294710433*x3-0.1358543740843192*x4+0.03668911137951071*x5);
      ch.F[97]=(-0.002868246193686554+0.00979803053682631*x-0.03456127641546138*x2+0.03441483523989725*x3-0.007427509325841953*x4-0.00210850167277718*x5);
      ch.F[98]=(0.001368517716503287-0.00616328670844153*x+0.0811980021841395*x2-0.1862323456516175*x3+0.1546505308397347*x4-0.04290213355892105*x5);
      ch.F[99]=(-0.002871667140060161-0.00825025391658457*x+0.04003409452300122*x2-0.06079587248994659*x3+0.03618410263920547*x4-0.007236388709495553*x5);
      ch.F[100]=(0.004088277269867671+0.006034387855027397*x-0.07991092992305502*x2+0.1849997128346744*x3-0.1550727516005535*x4+0.04341323398302099*x5);
      ch.F[101]=(-0.002421773572163585+0.006755512651017649*x-0.0454843112160964*x2+0.0870550686106944*x3-0.06499921386882105*x4+0.01667172223960862*x5);

}

