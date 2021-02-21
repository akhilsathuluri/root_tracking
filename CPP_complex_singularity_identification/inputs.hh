#ifndef INPUTS_H
#define INPUTS_H

#include<iostream>
#include<eigen3/Eigen/Dense>

using namespace Eigen;
// Add any input of solutions or leg lengths here and call them in the main function
/*!
An array of leg values corresponding to the path to be tracked in the task-space
*/
const double legvals_data[] = {1.4458208723874426,
   1.4060252316390276,
   1.3708899598114825,
   1.3403156514115564,
   1.3132546895462032,
   1.2881484000129075,
   1.263443138859881,
   1.2380553048623144,
   1.211745268649149,
   1.185388299791834,
   1.1610732176018594,
   1.1418722234161967,
   1.1311320840917536,
   1.1313661511986288,
   1.1432107771774895,
   1.1650585693451805,
   1.1935872752440149,
   1.2248021573719514,
   1.2550126767829055,
   1.2814337246761802,
   1.3024334715438375,
   1.3175791912107189,
   1.3275829446473018,
   1.3341408835158517,
   1.3395927442126638,
   1.3463695596368885,
   1.356351124710815,
   1.3704022458423257,
   1.38832306204365,
   1.4092141706965846,
   1.4320259156883506,
   1.4560198698916078,
   1.4809871670462658,
   1.5071977387443252,
   1.5351244027325475,
   1.5650237495550017,
   1.5965058230797666,
   1.6282735799507622,
   1.658188646424637,
   1.6836755810855188,
   1.7022879555614236,
   1.7121803954479813,
   1.7123255090072766,
   1.7024866163420196,
   1.683069475196287,
   1.6549821733884242,
   1.6195708145003058,
   1.5786204210694783,
   1.5343440888187254,
   1.489251214647897,
   1.4458208723874426,
   1.568679047566447,
   1.5698603026795732,
   1.5717537721590564,
   1.5731650379454716,
   1.5726066008225281,
   1.568597529116758,
   1.5598431815903437,
   1.5453140163790533,
   1.5242991496985907,
   1.4965057537674302,
   1.462227134595034,
   1.4225297323549373,
   1.3793398277837687,
   1.335295291530446,
   1.2933155479545482,
   1.2560119082445755,
   1.2251895676141775,
   1.2016615000138375,
   1.1854197232441477,
   1.176031384240781,
   1.1730577072601807,
   1.176326536784823,
   1.185958884090406,
   1.2021273711896059,
   1.2246271468186232,
   1.2524595032937231,
   1.2836758266949604,
   1.3156109115465828,
   1.3454062626863312,
   1.3705781982568874,
   1.3894317600928279,
   1.4012736122201563,
   1.4064860381329416,
   1.406526321840845,
   1.4038398024840395,
   1.4015872298272118,
   1.4030714768923358,
   1.4108853342135725,
   1.4260672795356202,
   1.4477183790961154,
   1.473342398908298,
   1.499720315554311,
   1.5238463443849297,
   1.5435557607830577,
   1.5577628253534497,
   1.5664206258378706,
   1.5703389080206165,
   1.5709291224985105,
   1.5698819456773074,
   1.5687764807078677,
   1.568679047566447,1.57865471573012,
   1.5772578276277807,
   1.5778550434070688,
   1.5804763044388688,
   1.5844641972877689,
   1.5888630480863817,
   1.5927561560814267,
   1.595446977387423,
   1.5964815651629782,
   1.5955644938018765,
   1.5924330814940793,
   1.5867498657306927,
   1.5780605576281015,
   1.5658375014045909,
   1.5495866105943508,
   1.528961666381093,1.50383320479944,
   1.4743000080842352,
   1.4406756438352373,
   1.40349862085779,
   1.3635976630869306,
   1.3222045854202875,
   1.2810559418775362,
   1.24237401197219,
   1.2086054101042907,
   1.1818819491644472,
   1.1633656544579651,
   1.1528157801790948,
   1.1486603317592439,
   1.148555947915099,
   1.1501524549925461,
   1.1517679103259142,
   1.1528440961679616,
   1.1541758422111095,
   1.1578953753598673,
   1.1671073725633816,
   1.185063064303903,
   1.2139773945843437,
   1.2539686839009128,
   1.3027423014606274,
   1.3562322726387697,
   1.4097848293457194,
   1.4592416713122303,
   1.5015573360610353,
   1.5349660686605489,
   1.5588967728694139,
   1.5738149545365978,
   1.5810640655554495,
   1.5826765171647856,
   1.5810785225308954,
   1.57865471573012,
   1.3393859531342835,
   1.368575630574579,
   1.3990402579634633,
   1.4288222696814576,
   1.4560204177960676,
   1.4792536849334956,
   1.4979004662831201,
   1.5121226533330083,
   1.5227379402845689,
   1.5309852696873718,
   1.538200251250867,
   1.5454304682575488,
   1.5530814556663117,
   1.5607387986983303,
   1.567277222847329,
   1.5712281052866268,
   1.5712373719657666,
   1.5664278762145638,
   1.5565814524800494,
   1.5421687593469104,
   1.5242967920395236,
   1.5046133362009289,
   1.4851492680829543,
   1.468043545906806,
   1.4551244591464778,
   1.4474256187466994,
   1.4448311580040185,
   1.4460444323070634,
   1.448908150047061,
   1.4509042076548107,
   1.449609006197701,
   1.4429911308797956,
   1.4295823954903513,
   1.4086219983903405,
   1.3802505848254434,
   1.3457450507661326,
   1.3076700025348424,
   1.269733252724591,
   1.2361677233858674,
   1.2106988517041266,
   1.1955019779514906,
   1.1906976421825644,
   1.1946459606372457,
   1.2048048595183738,
   1.2186701808472291,
   1.2344455213081487,
   1.2513377907280816,
   1.2695107554340237,
   1.2897514909301053,
   1.3129260105501703,
   1.3393859531342835,
   1.152479064455009,
   1.1595753765506795,
   1.1721049992320172,
   1.18887065543626,
   1.2082237126188393,
   1.2287137112670097,
   1.2496153967682926,
   1.271205099575883,
   1.294755760043515,1.32223440782205,
   1.3556932835577098,
   1.39644398567447,
   1.4442926386071655,
   1.4972251643060046,
   1.551765418849317,
   1.6038422568459851,
   1.6497250895566526,
   1.6866458219427636,
   1.7129988860704202,
   1.728242084930226,
   1.7326875702832472,
   1.7273157444390375,
   1.7136498271758644,
   1.6936489673279618,
   1.6695400137505436,
   1.643533013863014,
   1.6174536287772314,
   1.5924224975496022,
   1.568726664307539,
   1.5459313915737352,
   1.5231548106995787,
   1.4993786035995285,
   1.4737102502578532,
   1.445581392985175,
   1.4149016475530765,
   1.3821720077434474,
   1.348516891487555,
   1.3155607836244794,
   1.2851041180302696,
   1.2586598005680696,
   1.237030707426415,
   1.2201240883552567,
   1.2070750031070991,
   1.1965923286574918,
   1.187374323174641,1.17847893833186,
   1.1695966986532988,
   1.1611935747797215,
   1.154467481718555,
   1.1510521884154628,
   1.152479064455009,
   1.286875291086054,
   1.3058097061451694,
   1.3282124617412332,
   1.3504796245654223,
   1.3692233542469587,
   1.3817554335147053,
   1.386320639522505,
   1.3821931794665863,
   1.3697650289970196,
   1.3506710908047719,
   1.3278655903434193,
   1.3054421811294783,
   1.2879748728406808,
   1.27937368575113,
   1.2816708299306843,
   1.294416757562172,
   1.3150661385423303,
   1.3400693303583646,
   1.3660265964756986,
   1.3904549898083338,
   1.4120921979464427,
   1.430861208929243,1.44761804576645,
   1.463739689938029,
   1.4805963056788591,
   1.4990200036491172,
   1.5189598125080592,
   1.5394789171783798,
   1.559077833967764,1.57615099355277,
   1.5893554869769675,
   1.5977829431774786,
   1.6009558121152316,
   1.5987316031269319,
   1.5911944867441836,
   1.5785812388458356,
   1.561254321091324,
   1.5397076612593057,
   1.5145750608784607,
   1.4866145773003316,
   1.456665922936953,
   1.4256063286923566,
   1.3943406361440123,
   1.3638457627512215,
   1.335259281176652,1.30996952184303,
   1.289635087646386,
   1.2760442639628007,
   1.2707508478933789,
   1.2745312491735887,
   1.286875291086054};

   /*!
   Leg values corresponding to the path to be tracked in the task-space typecasted into a matrix
   */
Matrix<double, 51, 6> legvals(legvals_data);
// MatrixXcd Clegvals(51,6);
// Clegvals = legvals.cast<std::complex<double>>();

// Removing initsols and initphi from the main file
const double initsols_data[] = {
  0.21833808165665786,-0.7906213868033382,0.7240859710802774,-0.7154994439464201,-0.6452575011458942,0,
  -0.11224708177891103,0.7453128304748228,1.4299633878665543,0.8300450230008743,-0.2868399561208049,
  -0.24735531615022688,0.8761914454555555,1.3561766375736213,0.8054121712576314,0.7196393797208634,
  0.10661311496251118,-0.03391249606302923,0.20000000000000112,0,1.2799999999999998,0.19999999999999996,0,0,
  0.003056244705087313,-0.06688548576654774,0.45696277640159066,0.5469562635040984,-0.09469005512919006,
  0.0034901115067618653,0.3362762851850709,0.2868906143657115,-0.021027360477035077,0.024784412575612605,
  -0.3953833637520812,-0.37979521343323663,-0.5136951579345258,-0.5199956193725788,-0.7987229831575837,
  -0.3337590445754842,-0.9516994614991616,0,-1.885474722853147,-2.65366482193082,2.782283527303565,
  2.89251489256385,-2.2052227902999824,-1.742941092611359,0.6160601733585822,0.6974690501752805,
  0.4197629659432657,0.5377001757469011,0.07050254829194261,-0.10394753137509125,0.5244336967602403,
  0.26850949397918333,1.0307988965821262,0.6037494337274671,-0.18856114122522924,0,0.19200440495122106,
  0.08892583707114052,0.6828422266917203,1.0712524507882133,0.5587642531431387,0.32979660408779493,
  0.2758192040473777,0.26974569330836556,-0.13921668983029753,-0.3564442780964989,-1.0169800563674098,
  -0.6281098269592473,0.20000000000000112,0,-1.2799999999999998,-0.19999999999999996,0,0,3.138536408884706,
  -3.0747071678232456,2.6846298771882027,2.594636390085695,-3.046902598460603,3.1381025420830313,
  0.3362762851850709,0.2868906143657115,-0.021027360477035077,0.024784412575612605,-0.3953833637520812,
  -0.37979521343323663,-0.5136951579345258,-0.5199956193725788,0.7987229831575837,0.3337590445754842,
  0.9516994614991616,0,-1.2561179307366463,-0.48792783165897324,0.3593091262862285,0.24907776102594312,
  -0.936369863289811,-1.398651560978434,0.6160601733585822,0.6974690501752805,0.4197629659432657,
  0.5377001757469011,0.07050254829194261,-0.10394753137509125,0.21833808165665786,-0.7906213868033382,
  -0.7240859710802774,0.7154994439464201,0.6452575011458942,0,-3.0293455718108824,2.3962798231149707,
  1.7116292657232388,2.311547630588919,-2.8547526974689883,-2.8942373374395665,0.8761914454555555,
  1.3561766375736213,0.8054121712576314,0.7196393797208634,0.10661311496251118,-0.03391249606302923,
  0.5244336967602403,0.26850949397918333,-1.0307988965821262,-0.6037494337274671,0.18856114122522924,0,
  2.9495882486385723,3.0526668165186526,2.4587504268980727,2.07034020280158,2.5828284004466546,
  2.8117960495019982,0.2758192040473777,0.26974569330836556,-0.13921668983029753,-0.3564442780964989,
  -1.0169800563674098,-0.6281098269592473
};

Matrix<double, 8, 18>initsols(initsols_data);


const double phii_data[] = {
  0.2,0,1.28,0.2,0,0,
   0.003056244705087313,-0.06688548576654774,0.45696277640159066,0.5469562635040984,
   -0.09469005512919006,0.0034901115067618653,0.3362762851850709,0.2868906143657115,
   -0.021027360477035077,0.024784412575612605,-0.3953833637520812,-0.37979521343323663
};

Matrix<double, 18, 1>phii(phii_data);



const std::complex<double> phii_data2[] = {
  std::complex<double>(1.8700438348728161,0),std::complex<double>(-0.5235870043873335,0),std::complex<double>(0.,-0.926898973258719),std::complex<double>(0.,-0.10569096861540776),
     std::complex<double>(0.,0.7154784037290103),0,std::complex<double>(1.5707963267948966,1.5262959577414468),
     std::complex<double>(1.5707963267948966,0.6913993396133373),std::complex<double>(1.5707963267948966,-0.1997919002830573),
     std::complex<double>(1.5707963267948966,-0.5230490884925065),std::complex<double>(1.5707963267948966,0.37718441615032156),
     std::complex<double>(1.5707963267948966,1.3985097660833001),std::complex<double>(0.5889695511742916,-1.7974114210661284e-16),
     std::complex<double>(0.6273584382908439,-9.456869853007833e-17),std::complex<double>(0.3845316700324741,8.803289786100487e-17),
     std::complex<double>(0.5910427650086764,-9.488696919310518e-17),std::complex<double>(0.1527945222318842,-1.6367096645383268e-17),
     std::complex<double>(-0.07498723590723597,2.5098272394581245e-17)
};
Matrix<std::complex<double>, 18, 1>phii2(phii_data2);

const std::complex<double> initsols_data2[] = {
  std::complex<double>(1.8700438348728161,0),std::complex<double>(-0.5235870043873335,0),
  std::complex<double>(0.,-0.926898973258719),std::complex<double>(0.,-0.10569096861540776),
  std::complex<double>(0.,0.7154784037290103),std::complex<double>(0,0),
  std::complex<double>(1.5707963267948966,1.526295957741447),
  std::complex<double>(1.5707963267948966,0.6913993396133373),
  std::complex<double>(1.5707963267948966,-0.1997919002830573),
  std::complex<double>(1.5707963267948966,-0.5230490884925065),
  std::complex<double>(1.5707963267948966,0.37718441615032144),
  std::complex<double>(1.5707963267948966,1.3985097660833),
  std::complex<double>(0.5889695511742917,4.8919532403342845e-18),
  std::complex<double>(0.6273584382908439,-9.456869853007833e-17),
  std::complex<double>(0.38453167003247346,-1.2475304943139138e-17),
  std::complex<double>(0.5910427650086758,-6.21291993862106e-17),
  std::complex<double>(0.15279452223188406,-1.091982750937909e-16),
  std::complex<double>(-0.07498723590723586,3.3415727902921217e-17),
  std::complex<double>(1.8700438348728161,0),std::complex<double>(-0.5235870043873335,0),
  std::complex<double>(0.,0.926898973258719),std::complex<double>(0.,0.10569096861540776),
  std::complex<double>(0.,-0.7154784037290103),std::complex<double>(0,0),
  std::complex<double>(1.5707963267948966,-1.5262959577414468),
  std::complex<double>(1.5707963267948966,-0.6913993396133374),
  std::complex<double>(1.5707963267948966,0.19979190028305735),
  std::complex<double>(1.5707963267948966,0.5230490884925063),
  std::complex<double>(1.5707963267948966,-0.37718441615032183),
  std::complex<double>(1.5707963267948966,-1.3985097660833001),
  std::complex<double>(0.5889695511742917,4.8919532403342845e-18),
  std::complex<double>(0.6273584382908439,-9.456869853007833e-17),
  std::complex<double>(0.38453167003247346,-1.2475304943139138e-17),
  std::complex<double>(0.5910427650086758,-6.21291993862106e-17),
  std::complex<double>(0.15279452223188406,-1.091982750937909e-16),
  std::complex<double>(-0.07498723590723586,3.3415727902921217e-17),
  std::complex<double>(-1.0720337007448544,0),std::complex<double>(1.6971275035884408,0),
  std::complex<double>(0.,-0.9685245785328371),std::complex<double>(0.,-0.6638886705593215),
  std::complex<double>(0.,-0.3422565386596551),std::complex<double>(0,0),
  std::complex<double>(-1.5707963267948966,-0.7945901478802487),
  std::complex<double>(3.141592653589793,-0.8912894843397277),
  std::complex<double>(-3.141592653589793,-0.4613913273608144),
  std::complex<double>(-3.141592653589793,-0.9488972607958297),
  std::complex<double>(-1.5707963267948966,1.0375066349422688),
  std::complex<double>(-1.5707963267948966,0.8550510777135619),
  std::complex<double>(-0.8917219148146074,-1.779508141240245e-16),
  std::complex<double>(-1.5707963267948966,1.0539776567249441),
  std::complex<double>(-1.5707963267948966,1.2322351127771412),
  std::complex<double>(-1.5707963267948966,0.30446259360157596),
  std::complex<double>(-1.104806402560546,-1.2290999435359489e-16),
  std::complex<double>(-1.1334525087741785,-8.527491887348433e-17),
  std::complex<double>(-1.0720337007448544,0),std::complex<double>(1.6971275035884408,0),
  std::complex<double>(0.,0.9685245785328371),std::complex<double>(0.,0.6638886705593215),
  std::complex<double>(0.,0.3422565386596551),std::complex<double>(0,0),
  std::complex<double>(-1.5707963267948966,0.7945901478802487),
  std::complex<double>(0.,0.8912894843397278),
  std::complex<double>(9.776202109362632e-33,0.46139132736081423),
  std::complex<double>(0.,0.9488972607958296),
  std::complex<double>(-1.5707963267948966,-1.037506634942269),
  std::complex<double>(-1.5707963267948966,-0.855051077713562),
  std::complex<double>(-0.8917219148146074,-1.779508141240245e-16),
  std::complex<double>(-1.5707963267948966,1.0539776567249441),
  std::complex<double>(-1.5707963267948966,1.2322351127771412),
  std::complex<double>(-1.5707963267948966,0.30446259360157596),
  std::complex<double>(-1.104806402560546,-1.2290999435359489e-16),
  std::complex<double>(-1.1334525087741785,-8.527491887348433e-17),
  std::complex<double>(0.21833808165665786,0),std::complex<double>(-0.7906213868033382,0),std::complex<double>(0.7240859710802774,0),
  std::complex<double>(-0.7154994439464201,0),std::complex<double>(-0.6452575011458942,0),std::complex<double>(0,0),
  std::complex<double>(-0.11224708177891105,8.564716782405033e-18),
  std::complex<double>(0.7453128304748228,-7.710443659222416e-17),
  std::complex<double>(1.4299633878665543,-1.0123476658003623e-16),
  std::complex<double>(0.8300450230008743,-9.023128643113658e-17),
  std::complex<double>(-0.2868399561208049,-1.5752415902198853e-16),
  std::complex<double>(-0.24735531615022688,3.2662052072539e-18),
  std::complex<double>(0.8761914454555555,-1.513323722322559e-16),
  std::complex<double>(1.3561766375736213,8.564981963680112e-18),
  std::complex<double>(0.8054121712576314,-9.389105998133905e-18),
  std::complex<double>(0.7196393797208633,-1.9422490753309285e-16),
  std::complex<double>(0.10661311496251118,-2.0400690716215815e-17),
  std::complex<double>(-0.03391249606302923,-1.0051084355746774e-16),
  std::complex<double>(0.20000000000000112,0),std::complex<double>(0,0),std::complex<double>(1.2799999999999998,0),std::complex<double>(0.19999999999999996,0),std::complex<double>(0,0),std::complex<double>(0,0),
  std::complex<double>(0.003056244705087313,-8.472689376696798e-17),
  std::complex<double>(-0.06688548576654774,-8.940981208170962e-17),
  std::complex<double>(0.45696277640159066,-1.8850945593738954e-16),
  std::complex<double>(0.5469562635040984,-7.593834224960543e-17),
  std::complex<double>(-0.09469005512919007,-9.497861103054969e-17),
  std::complex<double>(0.0034901115067618653,-7.942309170426894e-17),
  std::complex<double>(0.3362762851850709,-3.591632290859678e-18),
  std::complex<double>(0.28689061436571145,-1.9452151002265927e-16),
  std::complex<double>(-0.021027360477035077,-1.3566684573936927e-16),
  std::complex<double>(0.02478441257561261,-1.7442214872157709e-16),
  std::complex<double>(-0.3953833637520812,-3.493768831479001e-17),
  std::complex<double>(-0.3797952134332366,-1.311572448543788e-16),
  std::complex<double>(-0.5136951579345258,0),std::complex<double>(-0.5199956193725788,0),std::complex<double>(-0.7987229831575837,0),
  std::complex<double>(-0.3337590445754842,0),std::complex<double>(-0.9516994614991616,0),std::complex<double>(0,0),
  std::complex<double>(-1.885474722853147,-1.0510305523517317e-16),
  std::complex<double>(-2.65366482193082,2.747262006499188e-17),
  std::complex<double>(2.782283527303565,1.860091889277488e-17),
  std::complex<double>(2.89251489256385,-1.5187686575730872e-16),
  std::complex<double>(-2.2052227902999824,5.3404586605507356e-17),
  std::complex<double>(-1.742941092611359,-1.4210208817509338e-16),
  std::complex<double>(0.6160601733585822,-1.1774174530588308e-16),
  std::complex<double>(0.6974690501752805,-2.1388868973451734e-16),
  std::complex<double>(0.4197629659432657,1.9431421639118244e-17),
  std::complex<double>(0.5377001757469011,-1.527409910691054e-16),
  std::complex<double>(0.07050254829194261,-6.627472181932422e-17),
  std::complex<double>(-0.10394753137509125,3.0597044027713674e-18),
  std::complex<double>(0.5244336967602403,0),std::complex<double>(0.26850949397918333,0),std::complex<double>(1.0307988965821262,0),
  std::complex<double>(0.6037494337274671,0),std::complex<double>(-0.18856114122522924,0),std::complex<double>(0,0),
  std::complex<double>(0.19200440495122106,-1.1303087763431644e-16),
  std::complex<double>(0.08892583707114052,-1.7414975994406066e-16),
  std::complex<double>(0.6828422266917203,-1.6291578756908741e-16),
  std::complex<double>(1.0712524507882133,2.178083104993402e-17),
  std::complex<double>(0.5587642531431386,-4.2184836469250337e-17),
  std::complex<double>(0.32979660408779493,-6.886554254330195e-17),
  std::complex<double>(0.2758192040473777,9.640119125577783e-18),
  std::complex<double>(0.26974569330836556,5.1519610043107026e-17),
  std::complex<double>(-0.13921668983029753,-2.7119620799532866e-17),
  std::complex<double>(-0.3564442780964989,-2.0790092391180416e-16),
  std::complex<double>(-1.0169800563674098,-1.8903285571184736e-16),
  std::complex<double>(-0.6281098269592473,-1.1081235020651129e-16),
  std::complex<double>(-0.7337396637589878,0),std::complex<double>(-2.1544767236034987,0),
  std::complex<double>(0.,-1.3894731322781502),std::complex<double>(0.,0.6338547336758182),
  std::complex<double>(0.,-0.42758251604766423),std::complex<double>(0,0),std::complex<double>(0.,0.4667881591550208),
  std::complex<double>(-1.5707963267948966,0.6723657450536872),
  std::complex<double>(-3.141592653589793,-0.1438830723200196),
  std::complex<double>(3.141592653589793,-0.3351950618101856),
  std::complex<double>(3.141592653589793,-0.6340224210835185),
  std::complex<double>(3.141592653589793,-1.039981928922122),
  std::complex<double>(1.5707963267948966,-0.23964948143486672),
  std::complex<double>(1.3880793950580876,-1.6589430343231204e-16),
  std::complex<double>(1.5707963267948966,-0.7214504540504836),
  std::complex<double>(1.5707963267948966,-1.6149681784690106),
  std::complex<double>(1.5707963267948966,-1.599193938876793),
  std::complex<double>(1.5707963267948966,-0.4520736617249368),
  std::complex<double>(-0.7337396637589878,0),std::complex<double>(-2.1544767236034987,0),
  std::complex<double>(0.,1.3894731322781502),std::complex<double>(0.,-0.6338547336758182),
  std::complex<double>(0.,0.42758251604766423),std::complex<double>(0,0),
  std::complex<double>(-3.141592653589793,-0.46678815915502114),
  std::complex<double>(-1.5707963267948966,-0.6723657450536874),
  std::complex<double>(1.4233381302634803e-32,0.14388307232001948),
  std::complex<double>(8.617153508413312e-33,0.3351950618101854),
  std::complex<double>(-1.1618308466671617e-32,0.6340224210835185),
  std::complex<double>(1.7436080963126474e-32,1.0399819289221217),
  std::complex<double>(1.5707963267948966,-0.23964948143486672),
  std::complex<double>(1.3880793950580876,-1.6589430343231204e-16),
  std::complex<double>(1.5707963267948966,-0.7214504540504836),
  std::complex<double>(1.5707963267948966,-1.6149681784690106),
  std::complex<double>(1.5707963267948966,-1.599193938876793),
  std::complex<double>(1.5707963267948966,-0.4520736617249368),
  std::complex<double>(0.20000000000000112,0),std::complex<double>(0,0),std::complex<double>(-1.2799999999999998,0),std::complex<double>(-0.19999999999999996,0),std::complex<double>(0,0),
  std::complex<double>(0,0),std::complex<double>(3.138536408884706,-8.472689376696798e-17),
  std::complex<double>(-3.0747071678232456,-8.940981208170962e-17),
  std::complex<double>(2.6846298771882027,-1.8850945593738954e-16),
  std::complex<double>(2.594636390085695,-7.593834224960543e-17),
  std::complex<double>(-3.046902598460603,-9.497861103054969e-17),
  std::complex<double>(3.1381025420830313,-7.942309170426894e-17),
  std::complex<double>(0.3362762851850709,-3.591632290859678e-18),
  std::complex<double>(0.28689061436571145,-1.9452151002265927e-16),
  std::complex<double>(-0.021027360477035077,-1.3566684573936927e-16),
  std::complex<double>(0.02478441257561261,-1.7442214872157709e-16),
  std::complex<double>(-0.3953833637520812,-3.493768831479001e-17),
  std::complex<double>(-0.3797952134332366,-1.311572448543788e-16),
  std::complex<double>(-0.5136951579345258,0),std::complex<double>(-0.5199956193725788,0),std::complex<double>(0.7987229831575837,0),
  std::complex<double>(0.3337590445754842,0),std::complex<double>(0.9516994614991616,0),std::complex<double>(0,0),
  std::complex<double>(-1.2561179307366463,-1.0510305523517317e-16),
  std::complex<double>(-0.4879278316589732,7.951949097344045e-17),
  std::complex<double>(0.3593091262862285,1.860091889277488e-17),
  std::complex<double>(0.24907776102594312,-1.5187686575730872e-16),
  std::complex<double>(-0.936369863289811,5.3404586605507356e-17),
  std::complex<double>(-1.398651560978434,-1.4210208817509338e-16),
  std::complex<double>(0.6160601733585822,-1.1774174530588308e-16),
  std::complex<double>(0.6974690501752805,-2.1388868973451734e-16),
  std::complex<double>(0.4197629659432657,1.9431421639118244e-17),
  std::complex<double>(0.5377001757469011,-1.527409910691054e-16),
  std::complex<double>(0.07050254829194261,-6.627472181932422e-17),
  std::complex<double>(-0.10394753137509125,3.0597044027713674e-18),
  std::complex<double>(0.21833808165665786,0),std::complex<double>(-0.7906213868033382,0),std::complex<double>(-0.7240859710802774,0),
  std::complex<double>(0.7154994439464201,0),std::complex<double>(0.6452575011458942,0),std::complex<double>(0,0),
  std::complex<double>(-3.0293455718108824,1.1673661063266579e-17),
  std::complex<double>(2.3962798231149702,-7.710443659222416e-17),
  std::complex<double>(1.7116292657232388,-1.0123476658003623e-16),
  std::complex<double>(2.311547630588919,-9.023128643113658e-17),
  std::complex<double>(-2.8547526974689883,-1.5752415902198853e-16),
  std::complex<double>(-2.8942373374395665,3.2662052072539e-18),
  std::complex<double>(0.8761914454555555,-1.513323722322559e-16),
  std::complex<double>(1.3561766375736213,8.564981963680112e-18),
  std::complex<double>(0.8054121712576314,-9.389105998133905e-18),
  std::complex<double>(0.7196393797208633,-1.9422490753309285e-16),
  std::complex<double>(0.10661311496251118,-2.0400690716215815e-17),
  std::complex<double>(-0.03391249606302923,-1.0051084355746774e-16),
  std::complex<double>(0.5244336967602403,0),std::complex<double>(0.26850949397918333,0),std::complex<double>(-1.0307988965821262,0),
  std::complex<double>(-0.6037494337274671,0),std::complex<double>(0.18856114122522924,0),std::complex<double>(0,0),
  std::complex<double>(2.9495882486385723,-1.1303087763431644e-16),
  std::complex<double>(3.0526668165186526,-1.7414975994406066e-16),
  std::complex<double>(2.4587504268980727,-1.6291578756908741e-16),
  std::complex<double>(2.0703402028015803,-3.1401649845148467e-17),
  std::complex<double>(2.5828284004466546,-4.2184836469250337e-17),
  std::complex<double>(2.8117960495019982,-6.886554254330195e-17),
  std::complex<double>(0.2758192040473777,9.640119125577783e-18),
  std::complex<double>(0.26974569330836556,5.1519610043107026e-17),
  std::complex<double>(-0.13921668983029753,-2.7119620799532866e-17),
  std::complex<double>(-0.3564442780964989,-2.0790092391180416e-16),
  std::complex<double>(-1.0169800563674098,-1.8903285571184736e-16),
  std::complex<double>(-0.6281098269592473,-1.1081235020651129e-16)
};
Matrix<std::complex<double>, 14, 18>initsols2(initsols_data2);

#endif