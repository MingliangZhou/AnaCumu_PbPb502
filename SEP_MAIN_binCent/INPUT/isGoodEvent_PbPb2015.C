int isGoodBCID_32(int bcid);
int isGoodBCID_33(int bcid);
int isGoodBCID_216(int bcid);
int isGoodBCID_400(int bcid);
int isGoodBCID_424(int bcid);
int isGoodBCID_492(int bcid);

int isGoodEvent_PbPb2015(int run, int lb, int bcid=1)
{
	switch(run)
	{
		case 287931:	if(lb>=143 && lb<=814)
						{
							return isGoodBCID_492(bcid);
						}
						else  return 0;

		case 287924:	if((lb>=75 && lb<=95) || (lb>=97 && lb<=367))
						{
							return isGoodBCID_492(bcid);
						}
						else  return 0;

		case 287866:	if(lb>=171 && lb<=793)
						{
							return isGoodBCID_492(bcid);
						}
						else  return 0;

		case 287843:	if((lb>=126 && lb<=184) || (lb>=186 && lb<=525))
						{
							return isGoodBCID_492(bcid);
						}
						else  return 0;

		case 287827:	if((lb>=51 && lb<=284) || (lb>=286 && lb<=458))
						{
							return isGoodBCID_492(bcid);
						}
						else  return 0;

		case 287728:	if(lb>=175 && lb<=577)
						{
							return isGoodBCID_492(bcid);
						}
						else  return 0;

		case 287706:	if(lb>=73 && lb<=440)
						{
							return isGoodBCID_492(bcid);
						}
						else  return 0;

		case 287632:	if(lb>=129 && lb<=497)
						{
							return isGoodBCID_424(bcid);
						}
						else  return 0;

		case 287594:	if(lb>=41 && lb<=194)
						{
							return isGoodBCID_424(bcid);
						}
						else  return 0;

		case 287560:	if((lb>=246 && lb<=269) || (lb>=272 && lb<=285) || (lb>=291 && lb<=295) || (lb>=304 && lb<=355))
						{
							return isGoodBCID_32(bcid);
						}
						else  return 0;

		case 287382:	if((lb>=135 && lb<=468) || (lb>=470 && lb<=477))
						{
							return isGoodBCID_424(bcid);
						}
						else  return 0;

		case 287380:	if(lb>=121 && lb<=124)
						{
							return isGoodBCID_424(bcid);
						}
						else  return 0;

		case 287378:	if((lb>=101 && lb<=270) || (lb>=272 && lb<=432))
						{
							return isGoodBCID_424(bcid);
						}
						else  return 0;

		case 287334:	if(lb>=333 && lb<=669)
						{
							return isGoodBCID_424(bcid);
						}
						else  return 0;

		case 287330:	if((lb>=99 && lb<=473) || (lb>=507 && lb<=538))
						{
							return isGoodBCID_424(bcid);
						}
						else  return 0;

		case 287321:	if(lb>=111 && lb<=170)
						{
							return isGoodBCID_424(bcid);
						}
						else  return 0;

		case 287281:	if(lb>=114 && lb<=584)
						{
							return isGoodBCID_424(bcid);
						}
						else  return 0;

		case 287270:	if(lb>=121 && lb<=576)
						{
							return isGoodBCID_424(bcid);
						}
						else  return 0;

		case 287259:	if(lb>=147 && lb<=613)
						{
							return isGoodBCID_400(bcid);
						}
						else  return 0;

		case 287224:	if((lb>=171 && lb<=195) || (lb>=197 && lb<=249))
						{
							return isGoodBCID_400(bcid);
						}
						else  return 0;

		case 287222:	if(lb>=261 && lb<=280)
						{
							return isGoodBCID_400(bcid);
						}
						else  return 0;

		case 287068:	if(lb>=500 && lb<=707)
						{
							return isGoodBCID_400(bcid);
						}
						else  return 0;

		case 287044:	if(lb>=251 && lb<=938)
						{
							return isGoodBCID_400(bcid);
						}
						else  return 0;

		case 287038:	if(lb>=143 && lb<=457)
						{
							return isGoodBCID_400(bcid);
						}
						else  return 0;

		case 286995:	return 0;

		case 286990:	if((lb>=183 && lb<=270) || (lb>=272 && lb<=314) || (lb>=316 && lb<=527))
						{
							return isGoodBCID_400(bcid);
						}
						else  return 0;

		case 286967:	return 0;

		case 286908:	if((lb>=144 && lb<=149) || (lb>=151 && lb<=180) || (lb>=182 && lb<=189) || (lb>=191 && lb<=199) || (lb>=201 && lb<=208) || (lb>=212 && lb<=215) || (lb>=217 && lb<=222) || (lb>=226 && lb<=240) || (lb>=242 && lb<=249) || (lb>=251 && lb<=434) || (lb>=436 && lb<=556))
						{
							return isGoodBCID_400(bcid);
						}
						else  return 0;

		case 286854:	if((lb>=433 && lb<=709) || (lb>=711 && lb<=915) || (lb>=917 && lb<=972))
						{
							return isGoodBCID_400(bcid);
						}
						else  return 0;

		case 286834:	if(lb>=153 && lb<=518)
						{
							return isGoodBCID_400(bcid);
						}
						else  return 0;		
		
		case 286767:	if(lb>=96 && lb<=395)
						{
							return isGoodBCID_216(bcid);
						}
						else  return 0;

		case 286748:	if((lb>=107 && lb<=155) || (lb>=161 && lb<=245) || (lb>=251 && lb<=368))
						{
							return isGoodBCID_216(bcid);
						}
						else  return 0;

		case 286717:	if(lb>=310 && lb<=477)
						{
							return isGoodBCID_33(bcid);
						}
						else  return 0;

		case 286711:	if(lb>=111 && lb<=235)
						{
							return isGoodBCID_33(bcid);
						}
						else  return 0;	
		
		case 286665:	if(lb>=156 && lb<=308)
						{
							return 1;
						}
						else  return 0;	
										
		default:		return -1;
	}
}

int isGoodBCID_32(int bcid)
{
	switch(bcid)
	{
		case 1:
		case 454:
		case 458:
		case 467:
		case 471:
		case 480:
		case 484:
		case 493:
		case 497:
		case 506:
		case 510:
		case 519:
		case 523:
		case 532:
		case 536:
		case 545:
		case 549:
		case 3119:
		case 3123:
		case 3132:
		case 3136:
		case 3145:
		case 3149:
		case 3158:
		case 3162:
		case 3171:
		case 3175:
		case 3184:
		case 3188:
		case 3197:
		case 3201:
		case 3210:
		case 3214:
					return 1;
		default:	return 0;
	}
}

int isGoodBCID_33(int bcid)
{
	switch(bcid)
	{
		case 1:
		case 454:
		case 458:
		case 467:
		case 471:
		case 480:
		case 484:
		case 493:
		case 497:
		case 506:
		case 510:
		case 519:
		case 523:
		case 532:
		case 536:
		case 545:
		case 549:
		case 3119:
		case 3123:
		case 3132:
		case 3136:
		case 3145:
		case 3149:
		case 3158:
		case 3162:
		case 3171:
		case 3175:
		case 3184:
		case 3188:
		case 3197:
		case 3201:
		case 3210:
		case 3214:
					return 1;
		default:	return 0;
	}
}

int isGoodBCID_216(int bcid)
{
	switch(bcid)
	{
		case 1:
		case 5:
		case 14:
		case 18:
		case 27:
		case 31:
		case 40:
		case 44:
		case 53:
		case 57:
		case 66:
		case 70:
		case 79:
		case 83:
		case 92:
		case 96:
		case 105:
		case 109:
		case 118:
		case 122:
		case 131:
		case 135:
		case 144:
		case 148:
		case 184:
		case 188:
		case 197:
		case 201:
		case 210:
		case 214:
		case 223:
		case 227:
		case 236:
		case 240:
		case 249:
		case 253:
		case 262:
		case 266:
		case 275:
		case 279:
		case 288:
		case 292:
		case 301:
		case 305:
		case 314:
		case 318:
		case 327:
		case 331:
		case 892:
		case 896:
		case 905:
		case 909:
		case 918:
		case 922:
		case 931:
		case 935:
		case 944:
		case 948:
		case 957:
		case 961:
		case 970:
		case 974:
		case 983:
		case 987:
		case 996:
		case 1000:
		case 1009:
		case 1013:
		case 1022:
		case 1026:
		case 1035:
		case 1039:
		case 1075:
		case 1079:
		case 1088:
		case 1092:
		case 1101:
		case 1105:
		case 1114:
		case 1118:
		case 1127:
		case 1131:
		case 1140:
		case 1144:
		case 1153:
		case 1157:
		case 1166:
		case 1170:
		case 1179:
		case 1183:
		case 1192:
		case 1196:
		case 1205:
		case 1209:
		case 1218:
		case 1222:
		case 1783:
		case 1787:
		case 1796:
		case 1800:
		case 1809:
		case 1813:
		case 1822:
		case 1826:
		case 1835:
		case 1839:
		case 1848:
		case 1852:
		case 1861:
		case 1865:
		case 1874:
		case 1878:
		case 1887:
		case 1891:
		case 1900:
		case 1904:
		case 1913:
		case 1917:
		case 1926:
		case 1930:
		case 1966:
		case 1970:
		case 1979:
		case 1983:
		case 1992:
		case 1996:
		case 2005:
		case 2009:
		case 2018:
		case 2022:
		case 2031:
		case 2035:
		case 2044:
		case 2048:
		case 2057:
		case 2061:
		case 2070:
		case 2074:
		case 2083:
		case 2087:
		case 2096:
		case 2100:
		case 2109:
		case 2113:
		case 2674:
		case 2678:
		case 2687:
		case 2691:
		case 2700:
		case 2704:
		case 2713:
		case 2717:
		case 2726:
		case 2730:
		case 2739:
		case 2743:
		case 2752:
		case 2756:
		case 2765:
		case 2769:
		case 2778:
		case 2782:
		case 2791:
		case 2795:
		case 2804:
		case 2808:
		case 2817:
		case 2821:
		case 2857:
		case 2861:
		case 2870:
		case 2874:
		case 2883:
		case 2887:
		case 2896:
		case 2900:
		case 2909:
		case 2913:
		case 2922:
		case 2926:
		case 2935:
		case 2939:
		case 2948:
		case 2952:
		case 2961:
		case 2965:
		case 2974:
		case 2978:
		case 2987:
		case 2991:
		case 3000:
		case 3004:
		case 3043:
		case 3047:
		case 3056:
		case 3060:
		case 3069:
		case 3073:
		case 3082:
		case 3086:
		case 3095:
		case 3099:
		case 3108:
		case 3112:
		case 3121:
		case 3125:
		case 3134:
		case 3138:
		case 3147:
		case 3151:
		case 3160:
		case 3164:
		case 3173:
		case 3177:
		case 3186:
		case 3190:
					return 1;
		default:	return 0;
	}
}

int isGoodBCID_400(int bcid)
{
	switch(bcid)
	{
		case 12:
		case 16:
		case 25:
		case 29:
		case 38:
		case 42:
		case 51:
		case 55:
		case 64:
		case 68:
		case 104:
		case 108:
		case 117:
		case 121:
		case 130:
		case 134:
		case 143:
		case 147:
		case 156:
		case 160:
		case 169:
		case 173:
		case 182:
		case 186:
		case 195:
		case 199:
		case 208:
		case 212:
		case 221:
		case 225:
		case 234:
		case 238:
		case 247:
		case 251:
		case 287:
		case 291:
		case 300:
		case 304:
		case 313:
		case 317:
		case 326:
		case 330:
		case 339:
		case 343:
		case 352:
		case 356:
		case 365:
		case 369:
		case 378:
		case 382:
		case 391:
		case 395:
		case 404:
		case 408:
		case 446:
		case 450:
		case 459:
		case 463:
		case 472:
		case 476:
		case 485:
		case 489:
		case 498:
		case 502:
		case 511:
		case 515:
		case 524:
		case 528:
		case 537:
		case 541:
		case 550:
		case 554:
		case 563:
		case 567:
		case 576:
		case 580:
		case 589:
		case 593:
		case 629:
		case 633:
		case 642:
		case 646:
		case 655:
		case 659:
		case 668:
		case 672:
		case 681:
		case 685:
		case 694:
		case 698:
		case 707:
		case 711:
		case 720:
		case 724:
		case 733:
		case 737:
		case 746:
		case 750:
		case 759:
		case 763:
		case 772:
		case 776:
		case 812:
		case 816:
		case 825:
		case 829:
		case 838:
		case 842:
		case 851:
		case 855:
		case 864:
		case 868:
		case 877:
		case 881:
		case 890:
		case 894:
		case 903:
		case 907:
		case 916:
		case 920:
		case 929:
		case 933:
		case 942:
		case 946:
		case 955:
		case 959:
		case 995:
		case 999:
		case 1008:
		case 1012:
		case 1021:
		case 1025:
		case 1034:
		case 1038:
		case 1047:
		case 1051:
		case 1060:
		case 1064:
		case 1073:
		case 1077:
		case 1086:
		case 1090:
		case 1099:
		case 1103:
		case 1112:
		case 1116:
		case 1125:
		case 1129:
		case 1138:
		case 1142:
		case 1178:
		case 1182:
		case 1191:
		case 1195:
		case 1204:
		case 1208:
		case 1217:
		case 1221:
		case 1230:
		case 1234:
		case 1243:
		case 1247:
		case 1256:
		case 1260:
		case 1269:
		case 1273:
		case 1282:
		case 1286:
		case 1520:
		case 1524:
		case 1533:
		case 1537:
		case 1546:
		case 1550:
		case 1559:
		case 1563:
		case 1572:
		case 1576:
		case 1585:
		case 1589:
		case 1598:
		case 1602:
		case 1611:
		case 1615:
		case 1624:
		case 1628:
		case 1637:
		case 1641:
		case 1650:
		case 1654:
		case 1663:
		case 1667:
		case 1703:
		case 1707:
		case 1716:
		case 1720:
		case 1729:
		case 1733:
		case 1742:
		case 1746:
		case 1755:
		case 1759:
		case 1768:
		case 1772:
		case 1781:
		case 1785:
		case 1794:
		case 1798:
		case 1807:
		case 1811:
		case 1820:
		case 1824:
		case 1833:
		case 1837:
		case 1846:
		case 1850:
		case 1886:
		case 1890:
		case 1899:
		case 1903:
		case 1912:
		case 1916:
		case 1925:
		case 1929:
		case 1938:
		case 1942:
		case 1951:
		case 1955:
		case 1964:
		case 1968:
		case 1977:
		case 1981:
		case 1990:
		case 1994:
		case 2003:
		case 2007:
		case 2016:
		case 2020:
		case 2029:
		case 2033:
		case 2069:
		case 2073:
		case 2082:
		case 2086:
		case 2095:
		case 2099:
		case 2108:
		case 2112:
		case 2121:
		case 2125:
		case 2134:
		case 2138:
		case 2147:
		case 2151:
		case 2160:
		case 2164:
		case 2173:
		case 2177:
		case 2186:
		case 2190:
		case 2228:
		case 2232:
		case 2241:
		case 2245:
		case 2254:
		case 2258:
		case 2267:
		case 2271:
		case 2280:
		case 2284:
		case 2293:
		case 2297:
		case 2306:
		case 2310:
		case 2319:
		case 2323:
		case 2332:
		case 2336:
		case 2345:
		case 2349:
		case 2358:
		case 2362:
		case 2371:
		case 2375:
		case 2411:
		case 2415:
		case 2424:
		case 2428:
		case 2437:
		case 2441:
		case 2450:
		case 2454:
		case 2463:
		case 2467:
		case 2476:
		case 2480:
		case 2489:
		case 2493:
		case 2502:
		case 2506:
		case 2515:
		case 2519:
		case 2528:
		case 2532:
		case 2541:
		case 2545:
		case 2554:
		case 2558:
		case 2594:
		case 2598:
		case 2607:
		case 2611:
		case 2620:
		case 2624:
		case 2633:
		case 2637:
		case 2646:
		case 2650:
		case 2659:
		case 2663:
		case 2672:
		case 2676:
		case 2685:
		case 2689:
		case 2698:
		case 2702:
		case 2711:
		case 2715:
		case 2724:
		case 2728:
		case 2737:
		case 2741:
		case 2777:
		case 2781:
		case 2790:
		case 2794:
		case 2803:
		case 2807:
		case 2816:
		case 2820:
		case 2829:
		case 2833:
		case 2842:
		case 2846:
		case 2855:
		case 2859:
		case 2868:
		case 2872:
		case 2881:
		case 2885:
		case 2894:
		case 2898:
		case 2907:
		case 2911:
		case 2920:
		case 2924:
		case 2960:
		case 2964:
		case 2973:
		case 2977:
		case 2986:
		case 2990:
		case 2999:
		case 3003:
		case 3012:
		case 3016:
		case 3025:
		case 3029:
		case 3038:
		case 3042:
		case 3051:
		case 3055:
		case 3064:
		case 3068:
		case 3077:
		case 3081:
		case 3119:
		case 3123:
		case 3132:
		case 3136:
		case 3145:
		case 3149:
		case 3158:
		case 3162:
		case 3171:
		case 3175:
		case 3184:
		case 3188:
		case 3197:
		case 3201:
		case 3210:
		case 3214:
		case 3223:
		case 3227:
		case 3236:
		case 3240:
		case 3249:
		case 3253:
		case 3262:
		case 3266:
					return 1;
		default:	return 0;
	}
}

int isGoodBCID_424(int bcid)
{
	switch(bcid)
	{
		case 1:
		case 18:
		case 22:
		case 29:
		case 33:
		case 40:
		case 44:
		case 51:
		case 55:
		case 62:
		case 66:
		case 73:
		case 77:
		case 84:
		case 88:
		case 124:
		case 128:
		case 135:
		case 139:
		case 146:
		case 150:
		case 157:
		case 161:
		case 168:
		case 172:
		case 179:
		case 183:
		case 190:
		case 194:
		case 201:
		case 205:
		case 212:
		case 216:
		case 223:
		case 227:
		case 234:
		case 238:
		case 245:
		case 249:
		case 285:
		case 289:
		case 296:
		case 300:
		case 307:
		case 311:
		case 318:
		case 322:
		case 329:
		case 333:
		case 340:
		case 344:
		case 351:
		case 355:
		case 362:
		case 366:
		case 373:
		case 377:
		case 384:
		case 388:
		case 395:
		case 399:
		case 607:
		case 611:
		case 618:
		case 622:
		case 629:
		case 633:
		case 640:
		case 644:
		case 651:
		case 655:
		case 693:
		case 697:
		case 704:
		case 708:
		case 715:
		case 719:
		case 726:
		case 730:
		case 737:
		case 741:
		case 748:
		case 752:
		case 759:
		case 763:
		case 770:
		case 774:
		case 781:
		case 785:
		case 792:
		case 796:
		case 803:
		case 807:
		case 814:
		case 818:
		case 854:
		case 858:
		case 865:
		case 869:
		case 876:
		case 880:
		case 887:
		case 891:
		case 898:
		case 902:
		case 909:
		case 913:
		case 920:
		case 924:
		case 931:
		case 935:
		case 942:
		case 946:
		case 953:
		case 957:
		case 964:
		case 968:
		case 975:
		case 979:
		case 1015:
		case 1019:
		case 1026:
		case 1030:
		case 1037:
		case 1041:
		case 1048:
		case 1052:
		case 1059:
		case 1063:
		case 1070:
		case 1074:
		case 1081:
		case 1085:
		case 1092:
		case 1096:
		case 1103:
		case 1107:
		case 1114:
		case 1118:
		case 1125:
		case 1129:
		case 1136:
		case 1140:
		case 1176:
		case 1180:
		case 1187:
		case 1191:
		case 1198:
		case 1202:
		case 1209:
		case 1213:
		case 1220:
		case 1224:
		case 1231:
		case 1235:
		case 1242:
		case 1246:
		case 1253:
		case 1257:
		case 1264:
		case 1268:
		case 1275:
		case 1279:
		case 1286:
		case 1290:
		case 1498:
		case 1502:
		case 1509:
		case 1513:
		case 1520:
		case 1524:
		case 1531:
		case 1535:
		case 1542:
		case 1546:
		case 1584:
		case 1588:
		case 1595:
		case 1599:
		case 1606:
		case 1610:
		case 1617:
		case 1621:
		case 1628:
		case 1632:
		case 1639:
		case 1643:
		case 1650:
		case 1654:
		case 1661:
		case 1665:
		case 1672:
		case 1676:
		case 1683:
		case 1687:
		case 1694:
		case 1698:
		case 1705:
		case 1709:
		case 1745:
		case 1749:
		case 1756:
		case 1760:
		case 1767:
		case 1771:
		case 1778:
		case 1782:
		case 1789:
		case 1793:
		case 1800:
		case 1804:
		case 1811:
		case 1815:
		case 1822:
		case 1826:
		case 1833:
		case 1837:
		case 1844:
		case 1848:
		case 1855:
		case 1859:
		case 1866:
		case 1870:
		case 1906:
		case 1910:
		case 1917:
		case 1921:
		case 1928:
		case 1932:
		case 1939:
		case 1943:
		case 1950:
		case 1954:
		case 1961:
		case 1965:
		case 1972:
		case 1976:
		case 1983:
		case 1987:
		case 1994:
		case 1998:
		case 2005:
		case 2009:
		case 2016:
		case 2020:
		case 2027:
		case 2031:
		case 2067:
		case 2071:
		case 2078:
		case 2082:
		case 2089:
		case 2093:
		case 2100:
		case 2104:
		case 2111:
		case 2115:
		case 2122:
		case 2126:
		case 2133:
		case 2137:
		case 2144:
		case 2148:
		case 2155:
		case 2159:
		case 2166:
		case 2170:
		case 2177:
		case 2181:
		case 2188:
		case 2192:
		case 2228:
		case 2232:
		case 2239:
		case 2243:
		case 2250:
		case 2254:
		case 2261:
		case 2265:
		case 2272:
		case 2276:
		case 2283:
		case 2287:
		case 2294:
		case 2298:
		case 2305:
		case 2309:
		case 2316:
		case 2320:
		case 2327:
		case 2331:
		case 2338:
		case 2342:
		case 2349:
		case 2353:
		case 2389:
		case 2393:
		case 2400:
		case 2404:
		case 2411:
		case 2415:
		case 2422:
		case 2426:
		case 2433:
		case 2437:
		case 2475:
		case 2479:
		case 2486:
		case 2490:
		case 2497:
		case 2501:
		case 2508:
		case 2512:
		case 2519:
		case 2523:
		case 2530:
		case 2534:
		case 2541:
		case 2545:
		case 2552:
		case 2556:
		case 2563:
		case 2567:
		case 2574:
		case 2578:
		case 2585:
		case 2589:
		case 2596:
		case 2600:
		case 2636:
		case 2640:
		case 2647:
		case 2651:
		case 2658:
		case 2662:
		case 2669:
		case 2673:
		case 2680:
		case 2684:
		case 2691:
		case 2695:
		case 2702:
		case 2706:
		case 2713:
		case 2717:
		case 2724:
		case 2728:
		case 2735:
		case 2739:
		case 2746:
		case 2750:
		case 2757:
		case 2761:
		case 2797:
		case 2801:
		case 2808:
		case 2812:
		case 2819:
		case 2823:
		case 2830:
		case 2834:
		case 2841:
		case 2845:
		case 2852:
		case 2856:
		case 2863:
		case 2867:
		case 2874:
		case 2878:
		case 2885:
		case 2889:
		case 2896:
		case 2900:
		case 2907:
		case 2911:
		case 2918:
		case 2922:
		case 2958:
		case 2962:
		case 2969:
		case 2973:
		case 2980:
		case 2984:
		case 2991:
		case 2995:
		case 3002:
		case 3006:
		case 3013:
		case 3017:
		case 3024:
		case 3028:
		case 3035:
		case 3039:
		case 3046:
		case 3050:
		case 3057:
		case 3061:
		case 3068:
		case 3072:
		case 3079:
		case 3083:
		case 3119:
		case 3123:
		case 3130:
		case 3134:
		case 3141:
		case 3145:
		case 3152:
		case 3156:
		case 3163:
		case 3167:
		case 3174:
		case 3178:
		case 3185:
		case 3189:
		case 3196:
		case 3200:
		case 3207:
		case 3211:
		case 3218:
		case 3222:
		case 3229:
		case 3233:
		case 3240:
		case 3244:
					return 1;
		default:	return 0;
	}
}

int isGoodBCID_492(int bcid)
{
	switch(bcid)
	{
		case 1:
		case 5:
		case 9:
		case 15:
		case 19:
		case 25:
		case 29:
		case 35:
		case 39:
		case 45:
		case 49:
		case 55:
		case 59:
		case 65:
		case 69:
		case 75:
		case 79:
		case 85:
		case 89:
		case 95:
		case 99:
		case 105:
		case 109:
		case 115:
		case 119:
		case 155:
		case 159:
		case 165:
		case 169:
		case 175:
		case 179:
		case 185:
		case 189:
		case 195:
		case 199:
		case 205:
		case 209:
		case 215:
		case 219:
		case 225:
		case 229:
		case 235:
		case 239:
		case 245:
		case 249:
		case 255:
		case 259:
		case 265:
		case 269:
		case 305:
		case 309:
		case 315:
		case 319:
		case 325:
		case 329:
		case 335:
		case 339:
		case 345:
		case 349:
		case 355:
		case 359:
		case 365:
		case 369:
		case 375:
		case 379:
		case 385:
		case 389:
		case 395:
		case 399:
		case 446:
		case 450:
		case 456:
		case 460:
		case 466:
		case 470:
		case 476:
		case 480:
		case 486:
		case 490:
		case 496:
		case 500:
		case 506:
		case 510:
		case 516:
		case 520:
		case 526:
		case 530:
		case 536:
		case 540:
		case 546:
		case 550:
		case 556:
		case 560:
		case 596:
		case 600:
		case 606:
		case 610:
		case 616:
		case 620:
		case 626:
		case 630:
		case 636:
		case 640:
		case 646:
		case 650:
		case 656:
		case 660:
		case 666:
		case 670:
		case 676:
		case 680:
		case 686:
		case 690:
		case 696:
		case 700:
		case 706:
		case 710:
		case 746:
		case 750:
		case 756:
		case 760:
		case 766:
		case 770:
		case 776:
		case 780:
		case 786:
		case 790:
		case 796:
		case 800:
		case 806:
		case 810:
		case 816:
		case 820:
		case 826:
		case 830:
		case 836:
		case 840:
		case 846:
		case 850:
		case 856:
		case 860:
		case 896:
		case 900:
		case 906:
		case 910:
		case 916:
		case 920:
		case 926:
		case 930:
		case 936:
		case 940:
		case 946:
		case 950:
		case 956:
		case 960:
		case 966:
		case 970:
		case 976:
		case 980:
		case 986:
		case 990:
		case 996:
		case 1000:
		case 1006:
		case 1010:
		case 1046:
		case 1050:
		case 1056:
		case 1060:
		case 1066:
		case 1070:
		case 1076:
		case 1080:
		case 1086:
		case 1090:
		case 1096:
		case 1100:
		case 1106:
		case 1110:
		case 1116:
		case 1120:
		case 1126:
		case 1130:
		case 1136:
		case 1140:
		case 1146:
		case 1150:
		case 1156:
		case 1160:
		case 1196:
		case 1200:
		case 1206:
		case 1210:
		case 1216:
		case 1220:
		case 1226:
		case 1230:
		case 1236:
		case 1240:
		case 1246:
		case 1250:
		case 1256:
		case 1260:
		case 1266:
		case 1270:
		case 1276:
		case 1280:
		case 1286:
		case 1290:
		case 1487:
		case 1491:
		case 1497:
		case 1501:
		case 1507:
		case 1511:
		case 1517:
		case 1521:
		case 1527:
		case 1531:
		case 1537:
		case 1541:
		case 1547:
		case 1551:
		case 1557:
		case 1561:
		case 1567:
		case 1571:
		case 1577:
		case 1581:
		case 1587:
		case 1591:
		case 1597:
		case 1601:
		case 1637:
		case 1641:
		case 1647:
		case 1651:
		case 1657:
		case 1661:
		case 1667:
		case 1671:
		case 1677:
		case 1681:
		case 1687:
		case 1691:
		case 1697:
		case 1701:
		case 1707:
		case 1711:
		case 1717:
		case 1721:
		case 1727:
		case 1731:
		case 1737:
		case 1741:
		case 1747:
		case 1751:
		case 1787:
		case 1791:
		case 1797:
		case 1801:
		case 1807:
		case 1811:
		case 1817:
		case 1821:
		case 1827:
		case 1831:
		case 1837:
		case 1841:
		case 1847:
		case 1851:
		case 1857:
		case 1861:
		case 1867:
		case 1871:
		case 1877:
		case 1881:
		case 1887:
		case 1891:
		case 1897:
		case 1901:
		case 1937:
		case 1941:
		case 1947:
		case 1951:
		case 1957:
		case 1961:
		case 1967:
		case 1971:
		case 1977:
		case 1981:
		case 1987:
		case 1991:
		case 1997:
		case 2001:
		case 2007:
		case 2011:
		case 2017:
		case 2021:
		case 2027:
		case 2031:
		case 2037:
		case 2041:
		case 2047:
		case 2051:
		case 2087:
		case 2091:
		case 2097:
		case 2101:
		case 2107:
		case 2111:
		case 2117:
		case 2121:
		case 2127:
		case 2131:
		case 2137:
		case 2141:
		case 2147:
		case 2151:
		case 2157:
		case 2161:
		case 2167:
		case 2171:
		case 2177:
		case 2181:
		case 2187:
		case 2191:
		case 2228:
		case 2232:
		case 2238:
		case 2242:
		case 2248:
		case 2252:
		case 2258:
		case 2262:
		case 2268:
		case 2272:
		case 2278:
		case 2282:
		case 2288:
		case 2292:
		case 2298:
		case 2302:
		case 2308:
		case 2312:
		case 2318:
		case 2322:
		case 2328:
		case 2332:
		case 2338:
		case 2342:
		case 2378:
		case 2382:
		case 2388:
		case 2392:
		case 2398:
		case 2402:
		case 2408:
		case 2412:
		case 2418:
		case 2422:
		case 2428:
		case 2432:
		case 2438:
		case 2442:
		case 2448:
		case 2452:
		case 2458:
		case 2462:
		case 2468:
		case 2472:
		case 2478:
		case 2482:
		case 2488:
		case 2492:
		case 2528:
		case 2532:
		case 2538:
		case 2542:
		case 2548:
		case 2552:
		case 2558:
		case 2562:
		case 2568:
		case 2572:
		case 2578:
		case 2582:
		case 2588:
		case 2592:
		case 2598:
		case 2602:
		case 2608:
		case 2612:
		case 2618:
		case 2622:
		case 2628:
		case 2632:
		case 2638:
		case 2642:
		case 2678:
		case 2682:
		case 2688:
		case 2692:
		case 2698:
		case 2702:
		case 2708:
		case 2712:
		case 2718:
		case 2722:
		case 2728:
		case 2732:
		case 2738:
		case 2742:
		case 2748:
		case 2752:
		case 2758:
		case 2762:
		case 2768:
		case 2772:
		case 2778:
		case 2782:
		case 2788:
		case 2792:
		case 2828:
		case 2832:
		case 2838:
		case 2842:
		case 2848:
		case 2852:
		case 2858:
		case 2862:
		case 2868:
		case 2872:
		case 2878:
		case 2882:
		case 2888:
		case 2892:
		case 2898:
		case 2902:
		case 2908:
		case 2912:
		case 2918:
		case 2922:
		case 2928:
		case 2932:
		case 2938:
		case 2942:
		case 2978:
		case 2982:
		case 2988:
		case 2992:
		case 2998:
		case 3002:
		case 3008:
		case 3012:
		case 3018:
		case 3022:
		case 3028:
		case 3032:
		case 3038:
		case 3042:
		case 3048:
		case 3052:
		case 3058:
		case 3062:
		case 3068:
		case 3072:
		case 3078:
		case 3082:
		case 3119:
		case 3123:
		case 3129:
		case 3133:
		case 3139:
		case 3143:
		case 3149:
		case 3153:
		case 3159:
		case 3163:
		case 3169:
		case 3173:
		case 3179:
		case 3183:
		case 3189:
		case 3193:
		case 3199:
		case 3203:
		case 3209:
		case 3213:
		case 3219:
		case 3223:
		case 3229:
		case 3233:
					return 1;
		default:	return 0;
	}
}
