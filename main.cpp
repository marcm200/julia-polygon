/*

	To construct polygons for the interior and exterior
	of trustworthily computed Julia sets and to answer
	the question, if a complex number belongs to one of those
	regions or cannot be determined
	
	Marc Meidlinger, August 2019
	
	GNU Public License GPL v3.0.
	No warranty.
	
	for the trustworthily computed Julia sets:
	see my repository: https://github.com/marcm200/julia-tsa-core
	
	Literature:
	
	The point-in-polygon test was performed by a simplified version of:
	A Simple and Correct Even-Odd Algorithm for the Point-in-Polygon Problem for Complex Polygons
	Conference Paper - February 2017
	Michael Galetzka, Patrick Glauner
	https://www.researchgate.net/publication/229158251

	RIGOROUS BOUNDS FOR POLYNOMIAL JULIA SETS
	Journal of Computational Dynamics 
	Luiz Henrique de Figueiredo1, Diego Nehab1
	Jorge Stolfi2 and Jo~ao Batista S. de Oliveira
	
	Images of Julia sets that you can trust
	LUIZ HENRIQUE DE FIGUEIREDO and DIEGO NEHAB
	JORGE STOLFI, JO ˜AO BATISTA OLIVEIRA

*/

#include <iostream>
#include "string.h"
#include "math.h"

typedef signed long long VLONG;
typedef unsigned char BYTE;

const BYTE COLORGRAY=0;
const BYTE COLORWHITE=0b01;
const BYTE COLORBLACK=0b10;
const BYTE COLORRED=4;
const BYTE COLORBLUE=5;
const BYTE EXTPOLCOL=COLORBLUE;
const BYTE COLORYELLOW=6;
const BYTE INTPOLCOL=COLORYELLOW;
const BYTE AKTIVCOL=16; 

const int MAXPOLYGONE=16384;
const int BORDERWIDTH=16;

enum { CMD_MAKEINT=1, CMD_MAKEEXT, CMD_QUALITY, CMD_ORACLE };
enum { PIP_ERROR=-1, PIP_UNKNOWN=0, PIP_INTERIOR, PIP_BOUNDARY, PIP_EXTERIOR };


// structs

struct RGB {
	BYTE R,G,B;
};

struct PlaneRect {
	double x0,x1,y0,y1;
};

struct Charmap {
	VLONG xlen,ylen; 
	VLONG memused;
	BYTE *cmp;
	RGB palette[256];
		
	Charmap();
	virtual ~Charmap();
		
	void fill(const BYTE);
	void copyFrom(Charmap&);
	void setlenxy(const int,const int);
	void saveAsBmp(const char*);
	int loadAsBmp(const char*);
	void setPaletteRGB(const int,const BYTE,const BYTE,const BYTE);
	inline void setPoint(const int,const int,const BYTE);
	inline BYTE getPoint(const int,const int);
	void lineVH(const int,const int,const int,const int,const BYTE);
	void fillrect(const int,const int,const int,const int,const BYTE);
};

struct PolygonPoint {
	int x,y;
};

struct Polygon {
	int pointcount;
	int memused;
	int useprepare;
	PolygonPoint *points;
	VLONG nenner; 
	double cx0,cx1,cy0,cy1;
	int xmin,xmax,ymin,ymax;
	int* yprepare;
		
	Polygon();
	virtual ~Polygon();
		
	void setlen(const int);
	int load(const char*);
	void save(const char*);
	void add(const int,const int);
	void trimColinearStart(void);
	int isColinearFree(void);
	int isDiagonalFree(void);
	void unPrepareY(void);
	void prepareY(const int);
};


// globals

Charmap inbild;
int granularity=5;
FILE *flog=NULL;
int RANGE0=-2,RANGE1=2;
int SCREENBREITE;
double skalaRangeProPixel;
int intpcount=0,extpcount=0;
Polygon *intp=NULL,*extp=NULL;
int LOWERBOUNDPOLYGONLENGTH=24;


// forward

#define LOGMSG(TT) \
{\
	printf(TT);\
	if (flog) fprintf(flog,TT);\
}

#define LOGMSG2(TT,TT2) \
{\
	printf(TT,TT2);\
	if (flog) fprintf(flog,TT,TT2);\
}

#define LOGMSG3(TT,TT2,TT3) \
{\
	printf(TT,TT2,TT3);\
	if (flog) fprintf(flog,TT,TT2,TT3);\
}

// functions corresponding to CMD_

int interiorPolygon(void);
int exteriorPolygon(void);
int jsoracle(const double,const double);
int qualitycontrol(void);

// constructing and testing functions

void oracle(const char*,const double,const double);
void oracleComplexNumber(const double,const double);
int point_in_polygonVH(Polygon&,const int,const int);
int qualitycontrol(Polygon& apg,const BYTE);
int buildPolygon(Charmap*,const char*);
Charmap* floodFillPattern(const int);
int preapreYOracle(const double);
void unPrepareYOracle(void);
int borderPresent(Charmap&);

// helper function
void loadAllPolygons(void);
void drawCrossing(Charmap*,const int,const int,const BYTE);
void drawAllPolygons(Charmap&);
void drawOnePolygon(Charmap&,Polygon&,const BYTE);
int inbildcoord(const double);


// small functions

inline int inbildcoord(const double w) {
	return (int)floor( (w - RANGE0) / skalaRangeProPixel );
}

inline int minimumI(const int a,const int b) {
	if (a < b) return a;
	return b;
}

inline int maximumI(const int a,const int b) {
	if (a > b) return a;
	return b;
}

inline void getMinMax(const int a,const int b,int& mi,int& ma) {
	if (a < b) { mi=a; ma=b; } else { mi=b; ma=a; }
}

char* chomp(char* s) {
	if (!s) return 0;
	for(int i=strlen(s);i>=0;i--) if (s[i]<32) s[i]=0; else break;
	return s;
}

char* upper(char* s) {
	if (!s) return 0;
	for(unsigned int i=0;i<strlen(s);i++) {
		if ((s[i]>='a')&&(s[i]<='z')) s[i]=s[i]-'a'+'A';
	}

	return s;
}


// struct Polygon

int Polygon::isDiagonalFree(void) {
	for(int i=1;i<pointcount;i++) {
		if (
			(points[i-1].x != points[i].x) &&
			(points[i-1].y != points[i].y) 
		) {
			return 0;
		}
	}
	
	return 1;
}

int Polygon::isColinearFree(void) {
	for(int i=2;i<pointcount;i++) {
		if (
			(
				(points[i-2].x == points[i].x) &&
				(points[i-1].x == points[i].x)
			) ||
			(
				(points[i-2].y == points[i].y) &&
				(points[i-1].y == points[i].y)
			) 
		) return 0;
	}
	
	// around the end check
	// [piintcount-1] is end point equal to [0] start point
	if (
		(
			(points[pointcount-2].x == points[1].x) &&
			(points[0].x == points[1].x)
		) ||
		(
			(points[pointcount-2].y == points[1].y) &&
			(points[0].y == points[1].y)
		) 
	) return 0;
	
	return 1;
}

void Polygon::trimColinearStart(void) {
	// identifies colinear segment around the end point
	
	while (pointcount >= 3) {
		if (
			(
				// x-kolinear
				(points[0].x == points[1].x) &&
				(points[0].x == points[pointcount-2].x)
			) ||
			(
				// y-kolinear
				(points[0].y == points[1].y) &&
				(points[0].y == points[pointcount-2].y)
			)
		) {
			// remove last point and set second to last
			// point as new start and end point
			pointcount--; 
			points[0].x=points[pointcount-1].x;
			points[1].x=points[pointcount-1].y;
		} else break;
	}
	
	if (pointcount <= 3) {
		LOGMSG("Possible error. colinear-trimming around end produced ver small untested polygon.\n");
	}
}

int Polygon::load(const char* afn) {
	FILE *f=fopen(afn,"rt");
	if (!f) return -1;
	char tmp[1024];
	int status=0;
	int punktnr=0;
	if (points) delete[] points;
	points=NULL;
	pointcount=0;
	int erster=1;
	
	nenner=0;
	fgets(tmp,1000,f); chomp(tmp);
	if (sscanf(tmp,"%I64d",&nenner) != 1) nenner=1 << 25;
	fgets(tmp,1000,f); chomp(tmp);
	if (sscanf(tmp,"%lf,%lf,%lf,%lf",&cx0,&cx1,&cy0,&cy1) != 4) {
		cx0=cy0=RANGE0;
		cx1=cy1=RANGE1;
	}
	fgets(tmp,1000,f); chomp(tmp);
	int a;
	if (sscanf(tmp,"%i",&a) != 1) {
		LOGMSG("ERROR. Polygon file not correct in point count.\n");
		exit(99);
	}
	setlen(a);
	pointcount=a;
	
	for(int i=0;i<pointcount;i++) {
		fgets(tmp,1000,f); chomp(tmp);
		int ax,ay;
		if (sscanf(tmp,"%i,%i",&ax,&ay) != 2) {
			LOGMSG2("ERROR. Polygon file not correct in point line %s.\n",tmp);
			exit(99);
		}
		points[i].x=ax;
		points[i].y=ay;
		if (i==0) {
			xmin=xmax=ax;
			ymin=ymax=ay;
		} else {
			if ( (ax-8) < xmin) xmin=ax-8;
			if ( (ax+8) > xmax) xmax=ax+8;
			if ( (ay-8) < ymin) ymin=ay-8;
			if ( (ay+8) > ymax) ymax=ay+8;
		}
	}

	fclose(f);

	return 1;
}

void Polygon::save(const char* afn) {
	FILE *f=fopen(afn,"wt");
	if (!f) return;
	
	fprintf(f,"%I64d\n",nenner);
	fprintf(f,"%i,%i,%i,%i\n",RANGE0,RANGE1,RANGE0,RANGE1);
	fprintf(f,"%i\n",pointcount);
	for(int i=0;i<pointcount;i++) {
		fprintf(f,"%i,%i\n",points[i].x,points[i].y);
	}
	
	fprintf(f,".\n");
	fclose(f);
}

Polygon::Polygon() {
	pointcount=0;
	points=NULL;
	memused=0;
	useprepare=0;
	yprepare=NULL;
}

Polygon::~Polygon() {
	if (points) delete[] points;
	if (yprepare) delete[] yprepare;
}

void Polygon::setlen(const int a) {
	if (points) delete[] points;
	memused=a;
	points=new PolygonPoint[memused];
	pointcount=0;
	useprepare=0;
	yprepare=new int[memused];
	// value not relevant before calling function setPrepareY
}

void Polygon::prepareY(const int ay) {
	useprepare=1;
	int li=-1;
	int BUFFER=2; // to account for rounding errors
	// "too many" intersections are considered valid

	for(int i=1;i<pointcount;i++) {
		if (
			(
				(points[i-1].y <= (ay+BUFFER)) &&
				(points[i].y >= (ay-BUFFER)) 
			) ||
			(
				(points[i-1].y >= (ay-BUFFER)) &&
				(points[i].y <= (ay+BUFFER)) 
			) 
		) {
			if (li >= 0) {
				yprepare[li]=i;
			} else {
				yprepare[0]=i; // first entry
			}
			li=i;
		}
	}
	
	useprepare=1;

	// jump to after the end of the polygon
	if (li>=0) {
		yprepare[li]=(pointcount+16);
	} else {
		yprepare[0]=(pointcount+16);
		// if no intersection occurs at all => nothing to test
	}
}

void Polygon::unPrepareY(void) {
	useprepare=0;
}

void Polygon::add(const int ax,const int ay) {
	// adds and checks directly for colinear segments
	
	if (pointcount > (memused-8)) {
		// memory is getting filled => allocate bigger part

		int mem2=memused+16384;
		PolygonPoint *p2=new PolygonPoint[mem2];
		for(int i=0;i<pointcount;i++) {
			p2[i].x=points[i].x;
			p2[i].y=points[i].y;
		}
		delete[] points;
		points=p2;
		memused=mem2;
	}
	
	if (pointcount == 0) {
		xmin=xmax=ax;
		ymin=ymax=ay;
	} else {
		if ( (ax-8) < xmin) xmin=ax-8;
		if ( (ax+8) > xmax) xmax=ax+8;
		if ( (ay-8) < ymin) ymin=ay-8;
		if ( (ay+8) > ymax) ymax=ay+8;
	}
	
	if (pointcount >= 2) {
		// ...A,B, new point
		if ( 
			(ax == points[pointcount-1].x) &&
			(ax == points[pointcount-2].x)
		) {
			// if colinear, then remove B
			points[pointcount-1].x=ax;
			points[pointcount-1].y=ay;
		} else
		if ( 
			(ay == points[pointcount-1].y) &&
			(ay == points[pointcount-2].y)
		) {
			points[pointcount-1].x=ax;
			points[pointcount-1].y=ay;
		} else {
			// point added
			points[pointcount].x=ax;
			points[pointcount].y=ay;
			pointcount++;
		}
	} else {
		points[pointcount].x=ax;
		points[pointcount].y=ay;
		pointcount++;
	}
}


// Charmap

void Charmap::fillrect(const int ax,const int ay,
	const int bx,const int by,const BYTE ff) {
	int x0=minimumI(ax,bx);
	int x1=maximumI(ax,bx);
	int y0=minimumI(ay,by);
	int y1=maximumI(ay,by);
	
	for(int y=y0;y<=y1;y++) for(int x=x0;x<=x1;x++) setPoint(x,y,ff);
	
}

void Charmap::copyFrom(Charmap& b) {
	if ( (xlen != b.xlen) || (ylen != b.ylen) ) return;
	if (!cmp) return;
	if (!b.cmp) return;
	
	for(int i=0;i<=255;i++) {
		palette[i].R=b.palette[i].R;
		palette[i].G=b.palette[i].G;
		palette[i].B=b.palette[i].B;
	}
	
	VLONG offset=0;
	for(int y=0;y<ylen;y++) for(int x=0;x<xlen;x++) {
		cmp[offset]=b.cmp[offset];
		offset++;
	}
}

void Charmap::setlenxy(const int ax,const int ay) {
	if ((xlen>0)&&(cmp)) delete[] cmp;
	
	xlen=ax;
	ylen=ay;
	memused=xlen*ylen;
	cmp=new BYTE[memused];
	if (!cmp) {
		LOGMSG("\nMemory error Charmap.\n");
		exit(99);
	}
}

Charmap::Charmap() {
	xlen=ylen=0;
	memused=0;
	cmp=NULL;
}

Charmap::~Charmap() {
	if ((xlen>0)&&(cmp)) delete[] cmp;
}

void Charmap::fill(const BYTE swert) {
	if (!cmp) return;
	for(VLONG i=0;i<memused;i++) cmp[i]=swert;
}

void Charmap::setPaletteRGB(const int pos,const BYTE ar,const BYTE ag,const BYTE ab) {
	if ((pos<0)||(pos>255)) return;
	palette[pos].R=ar;
	palette[pos].G=ag;
	palette[pos].B=ab;
}

void Charmap::setPoint(const int ax,const int ay,const BYTE awert) {
	VLONG pos=(VLONG)ay*xlen+ax;
	cmp[pos]=awert;
}

BYTE Charmap::getPoint(const int ax,const int ay) {
	VLONG pos=(VLONG)ay*xlen+ax;
	return cmp[pos];
}

void Charmap::lineVH(const int aax,const int aay,const int bbx,const int bby,const BYTE awert) {
	if (!cmp) return;
	
	int ax=aax;
	if (ax<0) ax=0;	else if (ax >= xlen) ax=xlen-1;
	int ay=aay;
	if (ay<0) ay=0;	else if (ay >= ylen) ay=ylen-1;
	int bx=bbx;
	if (bx<0) bx=0;	else if (bx >= xlen) bx=xlen-1;
	int by=bby;
	if (by<0) by=0;	else if (by >= ylen) by=ylen-1;
	
	if (ax == bx) {
		// vertical
		int y0,y1;
		if (ay < by) { y0=ay; y1=by; } else { y0=by; y1=ay; }
		for(int y=y0;y<=y1;y++) {
			VLONG pos=y*(VLONG)xlen+ax;
			cmp[pos]=awert;
		} // y
	} else if (ay == by) {
		// horizontal
		int x0,x1;
		if (ax < bx) { x0=ax; x1=bx; } else { x0=bx; x1=ax; }
		for(int x=x0;x<=x1;x++) {
			VLONG pos=ay*(VLONG)xlen+x;
			cmp[pos]=awert;
		} // y
	} 
	
	return;
}

unsigned char dez(const char c) {
	if ((c>='A')&&(c<='F')) { return c-'A'+10; }
	if ((c>='a')&&(c<='f')) { return c-'a'+10; }
	if ((c>='0')&&(c<='9')) { return c-'0'; }
	return 0;
}

void writehex(FILE *f,const char* s) {
	for(unsigned int i=0;i<strlen(s);i+=2) {
		unsigned char c = 16*dez(s[i]) + dez(s[i+1]);
		fwrite(&c,sizeof(c),1,f);
	}
}

void Charmap::saveAsBmp(const char* afn) {
	// 8 bit bitmap

	FILE *fbmp=fopen(afn,"wb");
	writehex(fbmp,"424D"); // BM
	int ybytes=ylen; 
	ybytes=(int)(4*ceil(ybytes*0.25));
	
	unsigned int off
		=		14 
			+	40 
			+	256*4
		;
	unsigned int filelen
			=	off
			+	(ybytes*xlen);
		;
			
	fwrite(&filelen,1,sizeof(filelen),fbmp);
	writehex(fbmp,"00000000"); 
	fwrite(&off,1,sizeof(off),fbmp); 
	writehex(fbmp,"28000000");
	
	unsigned int w = xlen;
	fwrite(&w,sizeof(w),1,fbmp);
	w = ylen;
	fwrite(&w,sizeof(w),1,fbmp);
	writehex(fbmp,"0100");
	writehex(fbmp,"0800");
	writehex(fbmp,"00000000");
	writehex(fbmp,"00000000");
	writehex(fbmp,"130B0000");
	writehex(fbmp,"130B0000");
	writehex(fbmp,"00010000");
	writehex(fbmp,"00000000");
	BYTE puffer[4];
	for(int i=0;i<256;i++) {
		puffer[0]=palette[i].B;
		puffer[1]=palette[i].G;
		puffer[2]=palette[i].R;
		puffer[3]=0;
		fwrite(puffer,4,sizeof(BYTE),fbmp);
	}
	
	fwrite(cmp,memused,sizeof(BYTE),fbmp);

	fclose(fbmp);	
}

int Charmap::loadAsBmp(const char* afn) {
	// must be 8-bit bitmap

	FILE *fbmp=fopen(afn,"rb");
	if (!fbmp) return 0;
	
	BYTE dummypuffer[4096];
#define DUMMYREAD(NR) \
	fread(dummypuffer,NR,sizeof(BYTE),fbmp);

	DUMMYREAD(2)
	int ybytes=ylen; 
	ybytes=(int)(4*ceil(ybytes*0.25));
	
	unsigned int off
		=		14 
			+	40 
			+	256*4 
		;
	unsigned int filelen
			=	off
			+	(ybytes*SCREENBREITE);
		;
			
	fread(&filelen,1,sizeof(filelen),fbmp);
	DUMMYREAD(4)
	fread(&off,1,sizeof(off),fbmp); 
	DUMMYREAD(4)
	
	unsigned int wx,wy;
	fread(&wx,sizeof(wx),1,fbmp);
	fread(&wy,sizeof(wy),1,fbmp);
	setlenxy(wx,wy);
	DUMMYREAD(2);
	unsigned short int bitsperpixel;
	fread(&bitsperpixel,sizeof(bitsperpixel),1,fbmp);
	if (bitsperpixel != 8) {
		LOGMSG("\n\nERROR. Image probably not 8-bit format.\n");
		exit(99);
	}
	
	DUMMYREAD(24)
	BYTE puffer[4];
	for(int i=0;i<256;i++) {
		fread(puffer,4,sizeof(BYTE),fbmp);
		palette[i].B=puffer[0];
		palette[i].G=puffer[1];
		palette[i].R=puffer[2];
	}
	
	fread(cmp,memused,sizeof(BYTE),fbmp);

	fclose(fbmp);
	
	return 1;	
}


// general functions

void setPaletteTo(Charmap& p) {
	for(int i=0;i<256;i++) {
		p.palette[i].R=0;
		p.palette[i].G=0;
		p.palette[i].B=0;
	}
	p.setPaletteRGB(COLORBLACK,0,0,0);
	p.setPaletteRGB(COLORWHITE,255,255,255);
	p.setPaletteRGB(COLORGRAY,127,127,127);
	p.setPaletteRGB(COLORRED,255,0,0);
	p.setPaletteRGB(COLORBLUE,0,0,255);
	p.setPaletteRGB(COLORYELLOW,255,255,0);
	p.setPaletteRGB(AKTIVCOL,0,255,127);
}

void adjustPalette(Charmap& md) {
	// take the RGB value of a pixel and set it
	// to the standard palette entries: COLORWHITE, etc...
	// so images from different sources can be used
		
	for(int y=0;y<md.ylen;y++) for(int x=0;x<md.xlen;x++) {
		int f=md.getPoint(x,y);
		if (
			(md.palette[f].R<20) &&
			(md.palette[f].G<20) &&
			(md.palette[f].B<20)
		) md.setPoint(x,y,COLORBLACK);
		else 
		if (
			(md.palette[f].R>230) &&
			(md.palette[f].G>230) &&
			(md.palette[f].B>230)
		) md.setPoint(x,y,COLORWHITE);
		else 
		if (
			(md.palette[f].R>50) &&
			(md.palette[f].G>50) &&
			(md.palette[f].B>50) &&
			(md.palette[f].R<200) &&
			(md.palette[f].G<200) &&
			(md.palette[f].B<200)
		) md.setPoint(x,y,COLORGRAY);
		else {
			LOGMSG("Error. Image contains invalid color.\n");
			exit(99);
		}
	}
	
	// now that the pixels have the correct palette
	// entry value, the palette is set accordingly
	setPaletteTo(md);
}

void drawCrossing(Charmap* inout,const int ax,const int ay,const BYTE af) {
	// red lines to indicte the region of error-causing pixels in an image
	inout->lineVH(0,ay-10,inout->xlen-1,ay-10,af);
	inout->lineVH(0,ay+10,inout->xlen-1,ay+10,af);
	inout->lineVH(ax-10,0,ax-10,inout->ylen-1,af);
	inout->lineVH(ax+10,0,ax+10,inout->ylen-1,af);
}

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//
// main construct function for interior polygons
//
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

int interiorPolygon(void) {
	// find interior regions and connect them
	// with space to the surrounding
	Charmap *blau=floodFillPattern(COLORBLACK);
	
	// find polygons and save them
	int erg=buildPolygon(blau,"int");
	
	delete blau;
	
	return erg;
}

// pattern-driven flood fill and connecting those
// within true interior or exterior as needed

Charmap* floodFillPattern(const int relf) {
	// look for GRANULARITYxGRANULARITY pattern of identical color
	// (black or white), mark the inner part leaving
	// one pixel at the circumference
	// then all interior or exterior has been marked
	// connect those kernels via appropriate
	// vertical or horizontal lines provided this line
	// does not touch any gray area
	// then those marked pixels which have as neighbour
	// at least one black or white (depending on what
	// type of polygon one is looking for) encompass the boundary

	Charmap *ptsa=new Charmap;
	ptsa->setlenxy(inbild.xlen,inbild.ylen);
	setPaletteTo(*ptsa);
	ptsa->copyFrom(inbild);
	
	if (granularity < 3) granularity=3;
	int D=granularity;
	printf("\nsearching for kernel points ...");
	for(int y=0;y<(inbild.ylen-D);y+=D) {
		for(int x=0;x<(inbild.xlen-D);x+=D) {
			if (ptsa->getPoint(x,y) != relf) continue;
				
			int gef=1;
			for(int dy=0;dy<D;dy++) {
				for(int dx=0;dx<D;dx++) {
					if (ptsa->getPoint(x+dx,y+dy) != relf) {
						gef=0;
						break;
					}
				} // dx
				if (gef <= 0) break;
			} // dy
				
			if (gef <= 0) continue;
				
			// leave border rectangle unchanged
			for(int y2=(y+1);y2<(y+D-1);y2++) {
				for(int x2=(x+1);x2<(x+D-1);x2++) {
					ptsa->setPoint(x2,y2,AKTIVCOL);
				}
			}
				
		} // x
	} // y
	
	// connect the patterns of AKTIVCOLOR
	int changed=1;
	printf("\nconnecting snippets ");
	while (changed>0) {
		printf(".");
		changed=0;
		for(int y=1;y<(ptsa->ylen-2);y++) {
			for(int x=1;x<(ptsa->xlen-2);x++) {
				if (ptsa->getPoint(x,y) != relf) continue;
			
				if (
					(ptsa->getPoint(x+1,y) == relf) &&
					(ptsa->getPoint(x-1,y) == AKTIVCOL) &&
					(ptsa->getPoint(x+2,y) == AKTIVCOL)
				) {
					// two horizontal black pixels with adjacent AKTIVCOL ones
					// line part above and below must be GRAY-free
					int korrekt=1;
					for(int dx=0;dx<4;dx++) {
						if (
							(ptsa->getPoint(x+dx,y-1) == COLORGRAY) ||
							(ptsa->getPoint(x+dx,y+1) == COLORGRAY)
						) {
							korrekt=0;
							break;
						}
					}
				
					if (korrekt>0) {
						// set the black pixels to AKTIVCOL as well to establish the connection
						ptsa->setPoint(x,y,AKTIVCOL);
						ptsa->setPoint(x+1,y,AKTIVCOL);
						changed=1;
					}
				} // gray active active gray
				else if (
					(ptsa->getPoint(x,y+1) == relf) &&
					(ptsa->getPoint(x,y+2) == AKTIVCOL) &&
					(ptsa->getPoint(x,y-1) == AKTIVCOL)
				) {
					// same in vertical connection
					int korrekt=1;
					for(int dy=0;dy<4;dy++) {
						if (
							(ptsa->getPoint(x-1,y+dy) == COLORGRAY) ||
							(ptsa->getPoint(x+1,y+dy) == COLORGRAY)
						) {
							korrekt=0;
							break;
						}
					}
				
					if (korrekt>0) {
						ptsa->setPoint(x,y,AKTIVCOL);
						ptsa->setPoint(x,y+1,AKTIVCOL);
						changed=1;
					}

				} 
			} // x
		} // y
	} // while changed
	
	// for exterior: AKTIVCOL can always to connected to the
	// screen border per requirement on
	// the image having BORDERWIDTH 
	if (relf == COLORWHITE) {
		ptsa->fillrect(0,0,ptsa->xlen-1,BORDERWIDTH-1,AKTIVCOL);
		ptsa->fillrect(0,ptsa->ylen-1-(BORDERWIDTH-1),ptsa->xlen-1,ptsa->ylen-1,AKTIVCOL);
		ptsa->fillrect(0,0,BORDERWIDTH-1,ptsa->ylen-1,AKTIVCOL);
		ptsa->fillrect(ptsa->xlen-1-(BORDERWIDTH-1),0,ptsa->xlen-1,ptsa->ylen-1,AKTIVCOL);
	}
	
	// boundary are those AKTIVCOL pixles with at least
	// one neighbour of RELF color
	printf("\nsearching for boundaries ...");
	
	for(int y=1;y<(ptsa->ylen-1);y++) {
		for(int x=1;x<(ptsa->xlen-1);x++) {
			if (ptsa->getPoint(x,y) != AKTIVCOL) continue;
			
			int gef=0;
			for(int dy=-1;dy<=1;dy++) {
				for(int dx=-1;dx<=1;dx++) {
					if (ptsa->getPoint(x+dx,y+dy) == relf) {
						gef=1;
						break;
					}
				}
				if (gef>0) break;
			}
			
			if (gef<=0) continue;
			
			ptsa->setPoint(x,y,COLORBLUE);
		} // x
	} // y
	
	// now the blue pixels constitute all the (to be found) polygons

	return ptsa;
}

int buildPolygon(Charmap* blau,const char* afnpref) {
	// looking repeatedly for blue pixels and following
	// them. By construction there should only ever be
	// two adjacent blue pixels to one blue pixel. So
	// the walking direction is determined (except for the first point in a polygon).
	
	VLONG NENNER=( (VLONG)1 << 25);
	int polanz=0;
	char tmp[1024];
	
	#define POLYGONADD(XX,YY) \
	{\
		int xp=(int)floor( ( (XX)*skalaRangeProPixel + RANGE0) * NENNER);\
		int yp=(int)floor( ( (YY)*skalaRangeProPixel + RANGE0) * NENNER);\
		p1->add(xp,yp);\
	}
	
	printf("\nsearching for polygons ");

	int changed=1;
	while (changed>0) {
		changed=0;
		printf(".");
		
		int startx=-1,starty=-1;
		
		// find an unused blue starting point of a polygon
		for(int y=0;y<blau->ylen;y++) {
			for(int x=0;x<blau->xlen;x++) {
				if (blau->getPoint(x,y) == COLORBLUE) {
					startx=x;
					starty=y;
					break;
				}
			}
			if (startx>=0) break;
		}
		
		if (startx < 0) {
			// no more starting points
			break; 
		}
		
		changed=1;
				
		// starting point has two unused blue neighbours
		// choose arbitrarily one to establish the direction
		int nx=-1,ny=-1,aktx=startx,akty=starty;
		if (blau->getPoint(aktx+1,akty) == COLORBLUE) {
			nx=aktx+1;
			ny=akty;
		} else
		if (blau->getPoint(aktx-1,akty) == COLORBLUE) {
			nx=aktx-1;
			ny=akty;
		} else
		if (blau->getPoint(aktx,akty-1) == COLORBLUE) {
			nx=aktx;
			ny=akty-1;
		} else
		if (blau->getPoint(aktx,akty+1) == COLORBLUE) {
			nx=aktx;
			ny=akty+1;
		} else {
			// one blue point surrounded by nothing
			// remove
			blau->setPoint(startx,starty,COLORYELLOW);
			continue; 
		}
		
		Polygon* p1=new Polygon;
		p1->setlen(blau->xlen << 4);
		p1->nenner=NENNER;
		p1->cx0=p1->cy0=RANGE0;
		p1->cx1=p1->cy1=RANGE1;
		blau->setPoint(nx,ny,COLORYELLOW);
		POLYGONADD(startx,starty);
		POLYGONADD(nx,ny);
		
		aktx=nx; akty=ny;
		int discard=0;
		
		// folow the current point aktx,akty to its
		// (should be) deterministic one unused blue neighbour
		// till the polygon is closed
		while (1) {
			if ( (aktx==startx) && (akty==starty) ) {
				// close polygon
				break;
			}
			
			nx=ny=-1;
			if (blau->getPoint(aktx+1,akty) == COLORBLUE) {
				nx=aktx+1;
				ny=akty;
			} else
			if (blau->getPoint(aktx-1,akty) == COLORBLUE) {
				nx=aktx-1;
				ny=akty;
			} else
			if (blau->getPoint(aktx,akty-1) == COLORBLUE) {
				nx=aktx;
				ny=akty-1;
			} else
			if (blau->getPoint(aktx,akty+1) == COLORBLUE) {
				nx=aktx;
				ny=akty+1;
			} else {
				// no further point to follow, but not closed
				// probably self-loop. Currently unhandled.
				LOGMSG("\n\nERROR. Polygon not closable. Probably self-loop.\n");
				drawCrossing(blau,aktx,akty,COLORRED);
				blau->saveAsBmp("_ERROR_not_closing.bmp");
				exit(99);
				discard=1;
				break;
			}
			
			// mark next point as visited
			blau->setPoint(nx,ny,COLORYELLOW);
			POLYGONADD(nx,ny)
			aktx=nx;
			akty=ny;
		} 
		
		if (discard<=0) {
			// valid poilygon. CHeck for colinearity
			// over the polygon's end.
			p1->trimColinearStart();
			if (p1->pointcount > LOWERBOUNDPOLYGONLENGTH) {
				sprintf(tmp,"%spoly%04i",afnpref,polanz);
				printf("possible polygon found with %i vertices: file %s\n",p1->pointcount,tmp);
				p1->save(tmp);
				polanz++;
			}
		} 
		
		delete p1;
	} // while
	
	return 1;
}

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//
// main routine to construct exterior polygons
// works the same way as for the interior ones except
// for a different color to look out for
//
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

int exteriorPolygon(void) {
	Charmap *blau=floodFillPattern(COLORWHITE);
	int erg=buildPolygon(blau,"ext");
	delete blau;
	
	return erg;
}

void calcSkala(void) {
	skalaRangeProPixel=(RANGE1-RANGE0); 
	skalaRangeProPixel /= SCREENBREITE;
}

void drawAllPolygons(Charmap& md) {
	for(int i=0;i<intpcount;i++) {
		drawOnePolygon(md,intp[i],COLORYELLOW);
	}
	for(int i=0;i<extpcount;i++) {
		drawOnePolygon(md,extp[i],COLORBLUE);
	}
}

int qcB2(
	Charmap& md,
	Polygon& apg,
	const BYTE relf,
	const BYTE apolcol
) {
	// check whether drawn-in polygons touch another
	// folfow each edge of a polygon and check
	// its neighbours
	
	#define COUNTNEIGHBOURS(XX,YY)\
	{\
		for(int dy=-1;dy<=1;dy++) {\
			for(int dx=-1;dx<=1;dx++) {\
				if ((dx==0)&&(dy==0)) continue;\
				BYTE f=md.getPoint(XX+dx,YY+dy);\
				if (f == relf) ctrrelf++;\
				else if (f== apolcol) ctrapolcol++;\
				else ctrother++;\
			}\
		}\
	}
	
	for(int i=1;i<apg.pointcount;i++) {
		// edge from [i-1] to [i]
		int xx0=inbildcoord( (double)apg.points[i-1].x / apg.nenner);
		int yy0=inbildcoord( (double)apg.points[i-1].y / apg.nenner);
		int xx1=inbildcoord( (double)apg.points[i].x / apg.nenner);
		int yy1=inbildcoord( (double)apg.points[i].y / apg.nenner);
		
		// end points: exactly TWO apolcol neighbous
		// other 6 have to be relf
		int ctrrelf,ctrapolcol,ctrother;
		ctrrelf=ctrapolcol=ctrother=0;
		COUNTNEIGHBOURS(xx0,yy0)
		if (
			(ctrapolcol != 2) || (ctrrelf != 6)
		) {
			LOGMSG("ERROR. Vertex wrong neighbours.\n");
			drawCrossing(&md,xx0,yy0,COLORRED);
			md.saveAsBmp("_ERROR_vertex.bmp");
			return 0;
		}
		ctrrelf=ctrapolcol=ctrother=0;
		COUNTNEIGHBOURS(xx1,yy1)
		if (
			(ctrapolcol != 2) || (ctrrelf != 6)
		) {
			LOGMSG("ERROR. Vertex wrong neighbours.\n");
			drawCrossing(&md,xx1,yy1,COLORRED);
			md.saveAsBmp("_ERROR_vertex.bmp");
			return 0;
		}

		if (xx0==xx1) {
			// vertical line
			int y0,y1;
			getMinMax(yy0,yy1,y0,y1);
			
			for(int y3=(y0+1);y3<=(y1-1);y3++) {
				// all points must be apolcol and left and right neighbour relf
				if (
					(md.getPoint(xx0-1,y3) == relf) &&
					(md.getPoint(xx0,y3) == apolcol) &&
					(md.getPoint(xx0+1,y3) == relf)
				) continue;
				else {
					LOGMSG("ERROR. Vertical line wrong.\n");
					drawCrossing(&md,xx0,y3,COLORRED);
					md.saveAsBmp("_ERROR_vertical.bmp");
					return 0;
				}
			}
		} else if (yy0==yy1) {
			// horizontal line
			int x0,x1;
			getMinMax(xx0,xx1,x0,x1);
			
			for(int x3=(x0+1);x3<=(x1-1);x3++) {
				// all points must be apolcol and upper and lower neighbour relf
				if (
					(md.getPoint(x3,yy0-1) == relf) &&
					(md.getPoint(x3,yy0) == apolcol) &&
					(md.getPoint(x3,yy0+1) == relf)
				) continue;
				else {
					LOGMSG("ERROR. Veritcal line wrong.\n");
					drawCrossing(&md,x3,yy0,COLORRED);
					md.saveAsBmp("_ERROR_vertical.bmp");
					return 0;
				}
			}
		}
	} // i
	return 1;
}

// polygons shall not cross each other
// and must have room to the surrounding border
int qcB(
	Charmap& md,
	Polygon& apg,
	const BYTE relf,
	const BYTE apolcol
) {
	
	// draw polygon's vertices and edges. They
	// are only allowed to hit pixels with color af /and LILA which the drawn points are set to temporarily)
	printf(".");
	int lx=-1,ly=-1;
	int encx0=md.xlen-1,encx1=0;
	int ency0=md.ylen-1,ency1=0;
	
	#define ENCL(XX,YY) \
	if ( (XX) < encx0 ) encx0=(XX);\
	if ( (XX) > encx1 ) encx1=(XX);\
	if ( (YY) < ency0 ) ency0=(YY);\
	if ( (YY) > ency1 ) ency1=(YY);
	
	// first pass: check whether polygon lies
	// completely in its region with space to
	// the surroundings
	for(int i=0;i<apg.pointcount;i++) {
		double d=apg.points[i].x; d /= apg.nenner;
		int xx=inbildcoord(d);
		d=apg.points[i].y; d /= apg.nenner;
		int yy=inbildcoord(d);
		if (lx>=0) {
			if (lx==xx) {
				// vertical line
				int y0,y1;
				getMinMax(ly,yy,y0,y1);
				for(int y3=y0;y3<=y1;y3++) {
					
					for(int dy=-1;dy<=1;dy++) {
						for(int dx=-1;dx<=1; dx++) {
							if (md.getPoint(xx+dx,y3+dy) != relf) {
								LOGMSG("ERROR. Polygon lies in wrong region.\n");
								drawCrossing(&md,xx,y3,COLORRED);
								md.saveAsBmp("_ERROR_wrong_region.bmp");
								return 0;
							}
						}
					}
				}
			} else if (ly==yy) {
				// horizontal line
				int x0,x1;
				getMinMax(lx,xx,x0,x1);
				for(int x3=x0;x3<=x1;x3++) {
					
					for(int dy=-1;dy<=1;dy++) {
						for(int dx=-1;dx<=1; dx++) {
							if (md.getPoint(x3+dx,yy+dy) != relf) {
								LOGMSG("ERROR. Polygon lies in wrong region.\n");
								drawCrossing(&md,x3,yy,COLORRED);
								md.saveAsBmp("_ERROR_wrong_region.bmp");
								return 0;
							}
						}
					}
				}
			} else {
				LOGMSG("ERROR. Diagonal.\n");
				drawCrossing(&md,lx,ly,COLORRED);
				drawCrossing(&md,xx,yy,COLORRED);
				md.saveAsBmp("_ERROR_diagonal.bmp");
				return 0;
			}
		} 
		
		lx=xx; 
		ly=yy;
	}
	
	printf(".");
	// now 2nd pass: draw the polygon in apolcol
	for(int i=0;i<apg.pointcount;i++) {
		double d=apg.points[i].x; d /= apg.nenner;
		int xx=inbildcoord(d);
		d=apg.points[i].y; d /= apg.nenner;
		int yy=inbildcoord(d);
		if (lx>=0) md.lineVH(lx,ly,xx,yy,apolcol);
		
		lx=xx;
		ly=yy;
	}
	
	return 1;
}

void drawOnePolygon(Charmap& md,Polygon& pol,const BYTE af) {
	double SKX=md.xlen; SKX /= (RANGE1-RANGE0);
	double SKY=md.ylen; SKY /= (RANGE1-RANGE0);
	
	int lx=-1,ly;
	
	for(int i=0;i<pol.pointcount;i++) {
		double d=pol.points[i].x; d /= pol.nenner;
		int sx=(int)floor( (d-RANGE0)*SKX );
		d=pol.points[i].y; d /= pol.nenner;
		int sy=(int)floor( (d-RANGE0)*SKY );
		
		if (lx>=0) md.lineVH(lx,ly,sx,sy,af);
		lx=sx;
		ly=sy;
	}
}

void unPrepareYOracle(void) {
	for(int i=0;i<intpcount;i++) intp[i].unPrepareY();
	for(int i=0;i<extpcount;i++) extp[i].unPrepareY();
}

void prepareYOracle(const double ay) {
	for(int i=0;i<intpcount;i++) {
		int py=(int)floor(ay*intp[i].nenner);
		intp[i].prepareY(py);
	}
	for(int i=0;i<extpcount;i++) {
		int py=(int)floor(ay*extp[i].nenner);
		extp[i].prepareY(py);
	}
}

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//
// main routine to answer whether a complex number
// lies within the set, or outside or cannot e determined
//
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

int jsoracle(const double ax,const double ay) {
	int mxy=2;
	// how many pixels in total with two layers
	// of neighbours
	int AREA=((mxy+mxy+1)*(mxy+mxy+1));
	
	for(int i=0;i<intpcount;i++) {
		int px=(int)floor(ax*intp[i].nenner);
		int py=(int)floor(ay*intp[i].nenner);
		int ic=0;
		for(int dy=-mxy;((ic>=0)&&(dy<=mxy));dy++) {
			for(int dx=-mxy;((ic>=0)&&(dx<=mxy));dx++) {
				if (point_in_polygonVH(intp[i],px+dx,py+dy) == PIP_INTERIOR) ic++;
				else ic=-1;
			}
		}
		
		// if in ONE interior polygon => return INTERIOR
		if (ic == AREA) return PIP_INTERIOR;
	} 
	
	int ergext=PIP_UNKNOWN;
	for(int i=0;i<extpcount;i++) if (extp[i].pointcount>0) {
		int ic=0;
		int px=(int)floor(ax*extp[i].nenner);
		int py=(int)floor(ay*extp[i].nenner);
		for(int dy=-mxy;((ic>=0)&&(dy<=mxy));dy++) {
			for(int dx=-mxy;((ic>=0)&&(dx<=mxy));dx++) {
				if (point_in_polygonVH(extp[i],px+dx,py+dy) == PIP_EXTERIOR) ic++;
				else ic=-1;
			}
		}
		
		if (ic != AREA) return PIP_UNKNOWN; else ergext=PIP_EXTERIOR;
	} 

	// if outside ALL etxerior polygons => exterior
	if (
		(ergext == PIP_EXTERIOR) && 
		(extpcount>0)
	) return PIP_EXTERIOR;

	return PIP_UNKNOWN;
}

void oracleComplexNumber(const double apx,const double apy) {
	// polygons must already be loaded
	if ( (intpcount<=0) && (extpcount<=0) ) {
		LOGMSG("\n\nERROR. No polygons loaded.\n");
		return;
	}
	LOGMSG3("point (%.20lg,%.20lg) ",apx,apy);
	int jserg=jsoracle(apx,apy);

	switch (jserg) {
		case PIP_INTERIOR: LOGMSG("definite INTERIOR\n"); break;
		case PIP_EXTERIOR: LOGMSG("definite EXTERIOR\n"); break;
		case PIP_UNKNOWN: LOGMSG("unknown\n"); break;
		default: LOGMSG2("\n\nERROR. jsoracle result %i\n",jserg); break;
	}
}

void oracle(const char* afn,const double apx,const double apy) {
	// load all polygons
	loadAllPolygons();

	if ((!afn) || (afn[0]<32)) {
		// one point
		oracleComplexNumber(apx,apy);
	} else {
		// check all points in the file
		FILE *f=fopen(afn,"rt");
		char tmp[1024];
		while (!feof(f)) {
			if (fgets(tmp,1000,f) <= 0) break;
			chomp(tmp);
			double px,py;
			if (sscanf(tmp,"%lf,%lf",&px,&py) == 2) {
				oracleComplexNumber(px,py);
			}
		}
		fclose(f);
	}
	
	delete[] intp;
	delete[] extp;
}

// test whether a (rational) point is inside or outsiude
// a closed polygon. The numbers ax,ay are to be interpreted as 
// rations with implicit denominator equal to variable "nenner" in the polygon
// this test is a simplified version of the even-odd
// test by Glauner et. al. as the polygons here only have vertical
// and horizontal parts
// "A Simple and Correct Even-Odd Algorithm for the Point-in-Polygon Problem for Complex Polygons"
// Conference Paper · February 2017
// Michael Galetzka, Patrick Glauner
// https://www.researchgate.net/publication/229158251

int point_in_polygonVH(
	Polygon& apg,
	const int ax,
	const int ay
) {
	// if outside the bounding rectangle of a polygon
	// the result is exterior
	
	if (
		(ax < apg.xmin) ||
		(ax > apg.xmax) ||
		(ay < apg.ymin) ||
		(ay > apg.ymax) 
	) return PIP_EXTERIOR;

	// if on the lines of the polygon (only horizontal
	// and vertical lines present by construction)
	// result is BOUNDARY
		
	int even=1;
	int i=0; // not 1, since incrementing is at start of while loop
	while (i<(apg.pointcount-1)) {
		if (apg.useprepare>0) {
			i=apg.yprepare[i];
		} else {
			i++;
		}
		if (i >= apg.pointcount) break;

		if (apg.points[i].x == apg.points[i-1].x) {
			// vertical line
			int miy=minimumI(apg.points[i].y,apg.points[i-1].y);
			int may=maximumI(apg.points[i].y,apg.points[i-1].y);
			
			if (
				(apg.points[i].x == ax) &&
				(miy <= ay) && (ay <= may)
			) return PIP_BOUNDARY;
			
			// does the horizontal ray intersect with this vertical segment
			if (
				(ax < apg.points[i].x) &&
				(miy <= ay) && 
				(ay <= may) 
			) {
				if (ay == apg.points[i].y) {
					int y0,y1,y2;
					y1=apg.points[i].y;
					if (i>0) y0=apg.points[i-1].y;
					else y0=apg.points[apg.pointcount-2].y;
					if (i<(apg.pointcount-1)) y2=apg.points[i+1].y;
					else y2=apg.points[1].y;
					
					if (
						( (y0 < y1) && (y1 < y2) ) ||
						( (y0 > y1) && (y1 > y2) )
					) {
						even=1-even;
					}
				} else 
				if (
					(miy < ay) && 
					(ay < may) 
				) even=1-even;
			}
		} else
		if (apg.points[i].y == apg.points[i-1].y) {
			// horizontal
			int minx=minimumI(apg.points[i].x,apg.points[i-1].x);
			int maxx=maximumI(apg.points[i].x,apg.points[i-1].x);
			
			if (
				(apg.points[i].y == ay) &&
				(minx <= ax) && (ax <= maxx)
			) return PIP_BOUNDARY;

			// are ray and segment colinear ?
			
			if (
				(ay == apg.points[i].y) &&
				(minx > ax)
			) {
				// does the ray enter or exit the polygon
				/*
				
					|  |	
					----	=> if ray touches the "----" line it does not enter the polygon
					
					|
					---	=> if ray colinearizes herw with "---" it changes its even od status
					  |		
				
				*/
				
				int y0,y1,y2;
				if (i > 1) {
					y0=apg.points[i-2].y;
					y1=apg.points[i].y;
				} else {
					if (i==1) {
						y0=apg.points[apg.pointcount-2].y;
						y1=apg.points[i].y;
					} else {
						y0=apg.points[apg.pointcount-2].y;
						y1=apg.points[0].y;
					}
				}
				
				if (i < (apg.pointcount-1)) {
					y2=apg.points[i+1].y;
				} else {
					y2=apg.points[1].y;
				}
				
				if (
					( (y0 < y1) && (y1 < y2) ) ||
					( (y0 > y1) && (y1 > y2) )
				) {
					even=1-even;
				}
				
			}
		} else {
			fprintf(flog,"\n\nERROR. Implementation. Diagonal #%i (%i,%i)->(%i,%i).\n",i-1,apg.points[i-1].x,apg.points[i-1].y,apg.points[i].x,apg.points[i].y);
			printf("\n\nERROR. Implementation. Diagonal #%i (%i,%i)->(%i,%i).\n",i-1,apg.points[i-1].x,apg.points[i-1].y,apg.points[i].x,apg.points[i].y);
			exit(99);
		}
	} // i
	
	if (even > 0) return PIP_EXTERIOR;
	
	return PIP_INTERIOR;
}

int qcA(Polygon& apg) {
	// closed
	
	if (
		(apg.points[0].x != apg.points[apg.pointcount-1].x) ||
		(apg.points[0].y != apg.points[apg.pointcount-1].y)
	) {
		LOGMSG("  ERROR: not closed\n");
		return 0;
	}

	// free of colinear segments
	// not necessary since it should be the case by construction

	if (apg.isColinearFree() <= 0) {
		LOGMSG("  ERROR. NOT free of colinear segments.\n");
		return 0;
	}
	
	// free of diagonals, not necessary either
	if (apg.isDiagonalFree() <= 0) {
		LOGMSG("  ERROR. NOT free of diagonal segments.\n");
		return 0;
	}
	
	return 1;
}

void loadAllPolygons(void) {
	if (extp) delete[] extp;
	if (intp) delete[] extp;
	
	extpcount=0;
	intpcount=0;
	extp=new Polygon[MAXPOLYGONE];
	intp=new Polygon[MAXPOLYGONE];
	
	int searche=1,searchi=1;
	char tmp[1024];
	
	while ( (searche>0) || (searchi>0) ) {
		if (searchi>0) {
			sprintf(tmp,"intpoly%04i",intpcount);
			if (intp[intpcount].load(tmp) <= 0) searchi=0; else {
				intpcount++;
			}
		}
		
		if (searche>0) {
			sprintf(tmp,"extpoly%04i",extpcount);
			if (extp[extpcount].load(tmp) <= 0) searche=0; else {
				extpcount++;
			}
		}
	}
}

int qualitycontrol(void) {
	int allvalid=1;
	Charmap small;
	int SMALLLEN=512;
	small.setlenxy(SMALLLEN,SMALLLEN);
	setPaletteTo(small);
	small.fill(COLORGRAY);
	
	// some non-power2 non-standard region
	double mm=0.5*(RANGE0+RANGE1);
	double br=0.783*(RANGE1-RANGE0);
	double sm0=mm-br;
	double sm1=mm+br;
	double smallskala=(double)(sm1-sm0) / SMALLLEN;
	
	loadAllPolygons();
	
	// Check A)
	int erg=1;
	LOGMSG("QC structure check: closed / colinear- and diagonal-free ... ");
	for(int i=0;i<intpcount;i++) if (qcA(intp[i]) <= 0) { erg=0; break; }
	if (erg>0) for(int i=0;i<extpcount;i++) if (qcA(extp[i]) <= 0) { erg=0; break; }
	if (erg <=0) {
		LOGMSG(" !! FAILED !!\n");
		return 0;
	}
	
	LOGMSG("\n  PASSED\n");

	// Check B
	LOGMSG("QC image check: positioning / spacing / cross- and touch-free ");
	for(int i=0;i<intpcount;i++) {
		if (qcB(inbild,intp[i],COLORBLACK,INTPOLCOL) <= 0) { 
			printf(" !! FAILED !!\n");
			return 0;
		}
	}
	
	for(int i=0;i<extpcount;i++) {
		if (qcB(inbild,extp[i],COLORWHITE,EXTPOLCOL) <= 0) { 
			printf(" !! FAILED !!\n");
			return 0;
		}
	}
	
	// now all polygons are drawn in their color
	// and they don't cross each other or touch
	// the boundary
	// do they touch one another ?
	
	
	printf(".");
	// go over all polygons again and follow their
	// edges. Check whether there are the right
	// count of polygon colored pixels and free ones
	// so that polygons do not touch each other
	for(int i=0;i<intpcount;i++) {
		if (qcB2(inbild,intp[i],COLORBLACK,INTPOLCOL) <= 0) {
			LOGMSG("FAILED.");
			return 0;
		}
	}
	for(int i=0;i<extpcount;i++) {
		if (qcB2(inbild,extp[i],COLORWHITE,EXTPOLCOL) <= 0) {
			LOGMSG("FAILED.");
			return 0;
		}
	}
		
	LOGMSG("\n  PASSED\n");

	// C-Test
	// bitmap-driven oracle test
	// every non-white pixel must lie on the exterior
	// of at least ONE exterior polygon
	// every non-black pixel must lie outside ALL
	// interior polygons
	
	int noch=1;
	int noch0;
	if (inbild.ylen <= 4096) noch0=inbild.ylen >> 3;
	else noch0=inbild.ylen >> 4;
	
	LOGMSG("QC oracle check: where do pixels lie with respect to polygon ");
	
	for(int y=0;y<inbild.ylen;y++) {
		double py=y*skalaRangeProPixel + RANGE0;
		if ( (--noch) <= 0) {
			printf("%I64lld ",inbild.ylen-y);
			noch=noch0;
		}
		
		for(int x=0;x<inbild.xlen;x++) {
			double px=x*skalaRangeProPixel + RANGE0;
			
			// exterior polygons
			if (inbild.getPoint(x,y) != COLORWHITE) {
				int tmp=intpcount;
				intpcount=0; 
				// temporary no interior polygons present
				// since only interested in functionality of exterior polygons here
				if (jsoracle(px,py) == PIP_EXTERIOR) {
					LOGMSG3("\n\nERROR. Exterior polygon tested wrong on image coordinates %i,%i\n",x,y);
					drawAllPolygons(inbild);
					drawCrossing(&inbild,x,y,COLORRED);
					inbild.saveAsBmp("_ERROR_quality.bmp");
					return 0;
				}
				intpcount=tmp;
			}
			
			// interier Polygons
			if (inbild.getPoint(x,y) != COLORBLACK) {
				int tmp=extpcount;
				extpcount=0; 
				// temporary no interior polygons
				if (jsoracle(px,py) == PIP_INTERIOR) {
					LOGMSG3("\n\nERROR. Interior polygon tested wrong on image coordinates %i,%i\n",x,y);
					drawAllPolygons(inbild);
					drawCrossing(&inbild,x,y,COLORRED);
					inbild.saveAsBmp("_ERROR_quality.bmp");
					return 0;
				}
				extpcount=tmp;
			}
		} // x
	} // y
	
	unPrepareYOracle();
	LOGMSG("\n  PASSED\n");
	LOGMSG("    i.e. no non-white pixel is judged as exterior\n");
	LOGMSG("    and  no non-black pixel is judged as interior\n");
	
	inbild.saveAsBmp("_FINAL_all_polygons.bmp");

	printf("\n\nadding to small image ...");
	double px,py;
	for(int y=0;y<SMALLLEN;y++) {
		py=(y+0.23)*smallskala + sm0;
		
		prepareYOracle(py);
		
		for(int x=0;x<SMALLLEN;x++) {
			if (small.getPoint(x,y) != COLORGRAY) continue;
			px=(x+0.23)*smallskala + sm0;
			int erg=jsoracle(px,py);
			
			switch (erg) {
				case PIP_EXTERIOR: small.setPoint(x,y,COLORWHITE); break;
				case PIP_INTERIOR: small.setPoint(x,y,COLORBLACK); break;
				case PIP_UNKNOWN: small.setPoint(x,y,COLORGRAY); break;
				default: LOGMSG("\n\nERROR. Small. jsOracle.\n"); return 0;
			}
		}
	}
	unPrepareYOracle();

	delete[] intp;
	delete[] extp;
	
	if (allvalid>0) {
		LOGMSG3("\n=========================================================\n\nVALID: Quality control: all consecutively numbered %i interior and %i exterior polygons passed the tests.\n\n=========================================================\n",intpcount,extpcount);
		small.saveAsBmp("_QC_passed_small_result.bmp");
		return 1;
	} else {
		LOGMSG("\nFAILURE: Quality control: set of polygons NOT USABLE.\n");
		return 0;
	}
	
	return allvalid;
}

int borderPresent(Charmap& md) {
	// image must have a white border. 
	int D=BORDERWIDTH;
	
	int border=1;
	// as image is quadratic, x/y can be interchanged
	for(int a=0;a<D;a++) {
		for(int b=0;b<md.xlen;b++) {
			if (
				(md.getPoint(a,b) != COLORWHITE) ||
				(md.getPoint(md.xlen-1-a,b) != COLORWHITE) ||
				(md.getPoint(b,a) != COLORWHITE) ||
				(md.getPoint(b,md.ylen-1-a) != COLORWHITE)
			) {
				border=0;
				break;
			}
		}
		
		if (border<=0) break;
	} // 
	
	return border;
}

// main entry
int main(int argc,char** argv) {
	flog=fopen("polygon.log","at");
	if (flog) fprintf(flog,"\n\n---------------\n");
	
	// standard values
	RANGE0=-2;
	RANGE1=2;
	granularity=5;
	int cmd=CMD_ORACLE;
	double px=0.0,py=0.0;
	char orakelfn[1024];
	orakelfn[0]=0;
	LOWERBOUNDPOLYGONLENGTH=24;
	
	// command line parameters
	// cmd=[makeint,makeext,quality,oracle]
	// range=a,b
	// granularity=n
	// minpollen=n
	// point=x,y or point=file
	
	for(int i=1;i<argc;i++) {
		upper(argv[i]);
		if (strstr(argv[i],"CMD=")==argv[i]) {
			if (strcmp(&argv[i][4],"MAKEINT")==0) cmd=CMD_MAKEINT;
			else if (strcmp(&argv[i][4],"MAKEEXT")==0) cmd=CMD_MAKEEXT;
			else if (strcmp(&argv[i][4],"ORACLE")==0) cmd=CMD_ORACLE;
			else if (strcmp(&argv[i][4],"QUALITY")==0) cmd=CMD_QUALITY;
		} else
		if (strstr(argv[i],"RANGE=")==argv[i]) {
			if (sscanf(&argv[i][6],"%i,%i",&RANGE0,&RANGE1) != 2) {
				RANGE0=-2;
				RANGE1=2;
			}
		} else
		if (strstr(argv[i],"POINT=")==argv[i]) {
			if (sscanf(&argv[i][6],"%lf,%lf",&px,&py) == 2) {
				orakelfn[0]=0;
			} else {
				if (strlen(argv[i])>1000) argv[i][1000]=0;
				strcpy(orakelfn,&argv[i][6]);
			}
		} else
		if (strstr(argv[i],"MINPOLLEN=")==argv[i]) {
			if (sscanf(&argv[i][10],"%i",&LOWERBOUNDPOLYGONLENGTH) != 1) {
				LOWERBOUNDPOLYGONLENGTH=24;
			}
		} else
		if (strstr(argv[i],"GRANULARITY=")==argv[i]) {
			if (sscanf(&argv[i][12],"%i",&granularity) != 1) {
				granularity=5;
			}
		} 
	} // i
	
	printf("loading image ...\n");
	if (inbild.loadAsBmp("_in.bmp") <= 0) {
		LOGMSG("\nERROR. Image _in.bmp not found.\n");
		exit(99);
	}
	adjustPalette(inbild);
	
	if (inbild.xlen != inbild.ylen) {
		LOGMSG("\nERROR. Only quadratic images feasible.\n");
		exit(99);
	}
	
	// image must have a white border of at least size 16
	if (borderPresent(inbild) <= 0) {
		LOGMSG("\nERROR. Image must have a white border.\n");
		exit(99);
	}
	
	SCREENBREITE=inbild.xlen;
	calcSkala();
	
	if (cmd==CMD_MAKEINT) interiorPolygon();
	else if (cmd==CMD_MAKEEXT) exteriorPolygon();
	else if (cmd==CMD_ORACLE) oracle(orakelfn,px,py);
	else if (cmd==CMD_QUALITY) qualitycontrol();
	
	if (flog) fclose(flog);

	return 0;
}

