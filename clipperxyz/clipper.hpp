/*******************************************************************************
*                                                                              *
* Author    :  Angus Johnson                                                   *
* Version   :  6.4.2                                                           *
* Date      :  27 February 2017                                                *
* Website   :  http://www.angusj.com                                           *
* Copyright :  Angus Johnson 2010-2017                                         *
*                                                                              *
* License:                                                                     *
* Use, modification & distribution is subject to Boost Software License Ver 1. *
* http://www.boost.org/LICENSE_1_0.txt                                         *
*                                                                              *
* Attributions:                                                                *
* The code in this library is an extension of Bala Vatti's clipping algorithm: *
* "A generic solution to polygon clipping"                                     *
* Communications of the ACM, Vol 35, Issue 7 (July 1992) pp 56-63.             *
* http://portal.acm.org/citation.cfm?id=129906                                 *
*                                                                              *
* Computer graphics and geometric modeling: implementation and algorithms      *
* By Max K. Agoston                                                            *
* Springer; 1 edition (January 4, 2005)                                        *
* http://books.google.com/books?q=vatti+clipping+agoston                       *
*                                                                              *
* See also:                                                                    *
* "Polygon Offsetting by Computing Winding Numbers"                            *
* Paper no. DETC2005-85513 pp. 565-575                                         *
* ASME 2005 International Design Engineering Technical Conferences             *
* and Computers and Information in Engineering Conference (IDETC/CIE2005)      *
* September 24-28, 2005 , Long Beach, California, USA                          *
* http://www.me.berkeley.edu/~mcmains/pubs/DAC05OffsetPolygon.pdf              *
*                                                                              *
*******************************************************************************/

#ifndef clipper_hpp_xyz
#define clipper_hpp_xyz

//use_int32: When enabled 32bit ints are used instead of 64bit ints. This
//improve performance but coordinate values are limited to the range +/- 46340
//#define use_int32

//use_lines: Enables line clipping. Adds a very minor cost to performance.
#define use_lines
  
//use_deprecated: Enables temporary support for the obsolete functions
//#define use_deprecated  

#include <vector>
#include <list>
#include <set>
#include <stdexcept>
#include <cstring>
#include <cstdlib>
#include <ostream>
#include <functional>
#include <queue>
#include <cmath>
namespace ClipperLibXYZ {

static double const pi = 3.141592653589793238;
static double const two_pi = pi * 2;
static double const def_arc_tolerance = 0.25;

static int const Unassigned = -1;  //edge not currently 'owning' a solution
static int const Skip = -2;        //edge that would otherwise close a path

#define HORIZONTAL (-1.0E+40)
#define TOLERANCE (1.0e-20)
#define NEAR_ZERO(val) (((val) > -TOLERANCE) && ((val) < TOLERANCE))

enum Direction { dRightToLeft, dLeftToRight };

enum ClipType { ctIntersection, ctUnion, ctDifference, ctXor };
enum PolyType { ptSubject, ptClip };
//By far the most widely used winding rules for polygon filling are
//EvenOdd & NonZero (GDI, GDI+, XLib, OpenGL, Cairo, AGG, Quartz, SVG, Gr32)
//Others rules include Positive, Negative and ABS_GTR_EQ_TWO (only in OpenGL)
//see http://glprogramming.com/red/chapter11.html
enum PolyFillType { pftEvenOdd, pftNonZero, pftPositive, pftNegative };

#ifdef use_int32
  typedef int cInt;
  static cInt const loRange = 0x7FFF;
  static cInt const hiRange = 0x7FFF;
#else
  typedef signed long long cInt;
  static cInt const loRange = 0x3FFFFFFF;
  static cInt const hiRange = 0x3FFFFFFFFFFFFFFFLL;
  typedef signed long long long64;     //used by Int128 class
  typedef unsigned long long ulong64;

#endif

struct IntPoint {
  cInt X;
  cInt Y;
  cInt Z;
  IntPoint(cInt x = 0, cInt y = 0, cInt z = 0): X(x), Y(y), Z(z) {};

  friend inline bool operator== (const IntPoint& a, const IntPoint& b)
  {
    return a.X == b.X && a.Y == b.Y;
  }
  friend inline bool operator!= (const IntPoint& a, const IntPoint& b)
  {
    return a.X != b.X  || a.Y != b.Y; 
  }
};

struct IntRect { cInt left; cInt top; cInt right; cInt bottom; };

//enums that are used internally ...
enum EdgeSide { esLeft = 1, esRight = 2 };

struct TEdge {
	IntPoint Bot;
	IntPoint Curr; //current (updated for every new scanbeam)
	IntPoint Top;
	double Dx;
	PolyType PolyTyp;
	EdgeSide Side; //side only refers to current side of solution poly
	int WindDelta; //1 or -1 depending on winding direction
	int WindCnt;
	int WindCnt2; //winding count of the opposite polytype
	int OutIdx;
	TEdge* Next;
	TEdge* Prev;
	TEdge* NextInLML;
	TEdge* NextInAEL;
	TEdge* PrevInAEL;
	TEdge* NextInSEL;
	TEdge* PrevInSEL;
};

struct IntersectNode {
	TEdge* Edge1;
	TEdge* Edge2;
	IntPoint        Pt;
};

struct LocalMinimum {
	cInt          Y;
	TEdge* LeftBound;
	TEdge* RightBound;
};

struct LocMinSorter
{
	inline bool operator()(const LocalMinimum& locMin1, const LocalMinimum& locMin2)
	{
		return locMin2.Y < locMin1.Y;
	}
};

TEdge* GetMaximaPair(TEdge* e);

#ifndef use_int32

//------------------------------------------------------------------------------
// Int128 class (enables safe math on signed 64bit integers)
// eg Int128 val1((long64)9223372036854775807); //ie 2^63 -1
//    Int128 val2((long64)9223372036854775807);
//    Int128 val3 = val1 * val2;
//    val3.AsString => "85070591730234615847396907784232501249" (8.5e+37)
//------------------------------------------------------------------------------

class Int128
{
public:
	ulong64 lo;
	long64 hi;

	Int128(long64 _lo = 0)
	{
		lo = (ulong64)_lo;
		if (_lo < 0)  hi = -1; else hi = 0;
	}


	Int128(const Int128& val) : lo(val.lo), hi(val.hi) {}

	Int128(const long64& _hi, const ulong64& _lo) : lo(_lo), hi(_hi) {}

	Int128& operator = (const long64& val)
	{
		lo = (ulong64)val;
		if (val < 0) hi = -1; else hi = 0;
		return *this;
	}
	Int128& operator = (const Int128& val) = default;

	bool operator == (const Int128& val) const
	{
		return (hi == val.hi && lo == val.lo);
	}

	bool operator != (const Int128& val) const
	{
		return !(*this == val);
	}

	bool operator > (const Int128& val) const
	{
		if (hi != val.hi)
			return hi > val.hi;
		else
			return lo > val.lo;
	}

	bool operator < (const Int128& val) const
	{
		if (hi != val.hi)
			return hi < val.hi;
		else
			return lo < val.lo;
	}

	bool operator >= (const Int128& val) const
	{
		return !(*this < val);
	}

	bool operator <= (const Int128& val) const
	{
		return !(*this > val);
	}

	Int128& operator += (const Int128& rhs)
	{
		hi += rhs.hi;
		lo += rhs.lo;
		if (lo < rhs.lo) hi++;
		return *this;
	}

	Int128 operator + (const Int128& rhs) const
	{
		Int128 result(*this);
		result += rhs;
		return result;
	}

	Int128& operator -= (const Int128& rhs)
	{
		*this += -rhs;
		return *this;
	}

	Int128 operator - (const Int128& rhs) const
	{
		Int128 result(*this);
		result -= rhs;
		return result;
	}

	Int128 operator-() const //unary negation
	{
		if (lo == 0)
			return Int128(-hi, 0);
		else
			return Int128(~hi, ~lo + 1);
	}

	operator double() const
	{
		const double shift64 = 18446744073709551616.0; //2^64
		if (hi < 0)
		{
			if (lo == 0) return (double)hi * shift64;
			else return -(double)(~lo + ~hi * shift64);
		}
		else
			return (double)(lo + hi * shift64);
	}

};

Int128 Int128Mul(long64 lhs, long64 rhs);

#endif
//------------------------------------------------------------------------------

typedef std::vector< IntPoint > Path;
typedef std::vector< Path > Paths;

inline Path& operator <<(Path& poly, const IntPoint& p) {poly.push_back(p); return poly;}
inline Paths& operator <<(Paths& polys, const Path& p) {polys.push_back(p); return polys;}

inline bool IsHorizontal(TEdge& e)
{
	return e.Dx == HORIZONTAL;
}
//------------------------------------------------------------------------------

inline double GetDx(const IntPoint pt1, const IntPoint pt2)
{
	return (pt1.Y == pt2.Y) ?
		HORIZONTAL : (double)(pt2.X - pt1.X) / (pt2.Y - pt1.Y);
}
//---------------------------------------------------------------------------

inline void SetDx(TEdge& e)
{
	cInt dy = (e.Top.Y - e.Bot.Y);
	if (dy == 0) e.Dx = HORIZONTAL;
	else e.Dx = (double)(e.Top.X - e.Bot.X) / dy;
}
//---------------------------------------------------------------------------

inline void SwapSides(TEdge& Edge1, TEdge& Edge2)
{
	EdgeSide Side = Edge1.Side;
	Edge1.Side = Edge2.Side;
	Edge2.Side = Side;
}
//------------------------------------------------------------------------------

inline void SwapPolyIndexes(TEdge& Edge1, TEdge& Edge2)
{
	int OutIdx = Edge1.OutIdx;
	Edge1.OutIdx = Edge2.OutIdx;
	Edge2.OutIdx = OutIdx;
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

inline cInt Round(double val)
{
	if ((val < 0)) return static_cast<cInt>(val - 0.5);
	else return static_cast<cInt>(val + 0.5);
}

inline cInt TopX(TEdge& edge, const cInt currentY)
{
	return (currentY == edge.Top.Y) ?
		edge.Top.X : edge.Bot.X + Round(edge.Dx * (currentY - edge.Bot.Y));
}

std::ostream& operator <<(std::ostream &s, const IntPoint &p);
std::ostream& operator <<(std::ostream &s, const Path &p);
std::ostream& operator <<(std::ostream &s, const Paths &p);

void save(const Paths& paths, const std::string& file);
void load(Paths& paths, const std::string& file);

struct DoublePoint
{
  double X;
  double Y;
  DoublePoint(double x = 0, double y = 0) : X(x), Y(y) {}
  DoublePoint(IntPoint ip) : X((double)ip.X), Y((double)ip.Y) {}
};
//------------------------------------------------------------------------------

typedef void (*ZFillCallback)(IntPoint& e1bot, IntPoint& e1top, IntPoint& e2bot, IntPoint& e2top, IntPoint& pt);

enum InitOptions {ioReverseSolution = 1, ioStrictlySimple = 2, ioPreserveCollinear = 4};
enum JoinType {jtSquare, jtRound, jtMiter};
enum EndType {etClosedPolygon, etClosedLine, etOpenButt, etOpenSquare, etOpenRound};

class PolyNode;
typedef std::vector< PolyNode* > PolyNodes;

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

inline cInt Abs(cInt val)
{
	return val < 0 ? -val : val;
}

class PolyNode 
{ 
public:
    PolyNode();
    virtual ~PolyNode(){};
    Path Contour;
    PolyNodes Childs;
    PolyNode* Parent;
    PolyNode* GetNext() const;
    bool IsHole() const;
    bool IsOpen() const;
    int ChildCount() const;
private:
    //PolyNode& operator =(PolyNode& other); 
    unsigned Index; //node index in Parent.Childs
    bool m_IsOpen;
    JoinType m_jointype;
    EndType m_endtype;
    PolyNode* GetNextSiblingUp() const;
    void AddChild(PolyNode& child);
    friend class Clipper; //to access Index
    friend class ClipperOffset; 
	friend class ClipperOffsetEx;
	friend class ClipperEx;
};

struct OutPt {
	int       Idx;
	IntPoint  Pt;
	OutPt* Next;
	OutPt* Prev;
};

//OutRec: contains a path in the clipping solution. Edges in the AEL will
//carry a pointer to an OutRec when they are part of the clipping solution.
struct OutRec {
	int       Idx;
	bool      IsHole;
	bool      IsOpen;
	OutRec* FirstLeft;  //see comments in clipper.pas
	PolyNode* PolyNd;
	OutPt* Pts;
	OutPt* BottomPt;
};



struct Join {
	OutPt* OutPt1;
	OutPt* OutPt2;
	IntPoint  OffPt;
};

OutRec* GetLowermostRec(OutRec* outRec1, OutRec* outRec2);
bool OutRec1RightOfOutRec2(OutRec* outRec1, OutRec* outRec2);
TEdge* GetMaximaPairEx(TEdge* e);

void SwapIntersectNodes(IntersectNode& int1, IntersectNode& int2);
bool GetOverlap(const cInt a1, const cInt a2, const cInt b1, const cInt b2,
	cInt& Left, cInt& Right);

bool JoinHorz(OutPt* op1, OutPt* op1b, OutPt* op2, OutPt* op2b,
	const IntPoint Pt, bool DiscardLeft);

TEdge* GetNextInAEL(TEdge* e, Direction dir);
//------------------------------------------------------------------------------

void GetHorzDirection(TEdge& HorzEdge, Direction& Dir, cInt& Left, cInt& Right);

bool IntersectListSort(IntersectNode* node1, IntersectNode* node2);
int PointCount(OutPt* Pts);

OutPt* DupOutPt(OutPt* outPt, bool InsertAfter);
//------------------------------------------------------------------------------

class PolyTree: public PolyNode
{ 
public:
    ~PolyTree(){ Clear(); };
    PolyNode* GetFirst() const;
    void Clear();
    int Total() const;
private:
  //PolyTree& operator =(PolyTree& other);
  PolyNodes AllNodes;
    friend class Clipper; //to access AllNodes
	friend class ClipperEx;
};

bool Orientation(const Path &poly);
double Area(const Path &poly);
int PointInPolygon(const IntPoint &pt, const Path &path);

void SimplifyPolygon(const Path &in_poly, Paths &out_polys, PolyFillType fillType = pftEvenOdd);
void SimplifyPolygons(const Paths &in_polys, Paths &out_polys, PolyFillType fillType = pftEvenOdd);
void SimplifyPolygons(Paths &polys, PolyFillType fillType = pftEvenOdd);

void CleanPolygon(const Path& in_poly, Path& out_poly, double distance = 1.415);
void CleanPolygon(Path& poly, double distance = 1.415);
void CleanPolygons(const Paths& in_polys, Paths& out_polys, double distance = 1.415);
void CleanPolygons(Paths& polys, double distance = 1.415);

void MinkowskiSum(const Path& pattern, const Path& path, Paths& solution, bool pathIsClosed);
void MinkowskiSum(const Path& pattern, const Paths& paths, Paths& solution, bool pathIsClosed);
void MinkowskiDiff(const Path& poly1, const Path& poly2, Paths& solution);

void PolyTreeToPaths(const PolyTree& polytree, Paths& paths);
void ClosedPathsFromPolyTree(const PolyTree& polytree, Paths& paths);
void OpenPathsFromPolyTree(PolyTree& polytree, Paths& paths);

void ReversePath(Path& p);
void ReversePaths(Paths& p);

bool Orientation(const Path& poly);
//------------------------------------------------------------------------------

double Area(const Path& poly);
//------------------------------------------------------------------------------

double Area(const OutPt* op);
//------------------------------------------------------------------------------

double Area(const OutRec& outRec);
//------------------------------------------------------------------------------

bool PointIsVertex(const IntPoint& Pt, OutPt* pp);
//------------------------------------------------------------------------------

//See "The Point in Polygon Problem for Arbitrary Polygons" by Hormann & Agathos
//http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.88.5498&rep=rep1&type=pdf
int PointInPolygon(const IntPoint& pt, const Path& path);
//------------------------------------------------------------------------------

int PointInPolygon(const IntPoint& pt, OutPt* op);
//------------------------------------------------------------------------------

bool Poly2ContainsPoly1(OutPt* OutPt1, OutPt* OutPt2);
//----------------------------------------------------------------------

bool SlopesEqual(const TEdge& e1, const TEdge& e2, bool UseFullInt64Range);
//------------------------------------------------------------------------------

bool SlopesEqual(const IntPoint pt1, const IntPoint pt2,
	const IntPoint pt3, bool UseFullInt64Range);
//------------------------------------------------------------------------------

bool SlopesEqual(const IntPoint pt1, const IntPoint pt2,
	const IntPoint pt3, const IntPoint pt4, bool UseFullInt64Range);
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------

void IntersectPoint(TEdge& Edge1, TEdge& Edge2, IntPoint& ip);
//------------------------------------------------------------------------------

void ReversePolyPtLinks(OutPt* pp);
//------------------------------------------------------------------------------

void DisposeOutPts(OutPt*& pp);
//------------------------------------------------------------------------------

inline void InitEdge(TEdge* e, TEdge* eNext, TEdge* ePrev, const IntPoint& Pt)
{
	std::memset(e, 0, sizeof(TEdge));
	e->Next = eNext;
	e->Prev = ePrev;
	e->Curr = Pt;
	e->OutIdx = Unassigned;
}
//------------------------------------------------------------------------------

void InitEdge2(TEdge& e, PolyType Pt);
//------------------------------------------------------------------------------

TEdge* RemoveEdge(TEdge* e);
//------------------------------------------------------------------------------

inline void ReverseHorizontal(TEdge& e)
{
	//swap horizontal edges' Top and Bottom x's so they follow the natural
	//progression of the bounds - ie so their xbots will align with the
	//adjoining lower edge. [Helpful in the ProcessHorizontal() method.]
	std::swap(e.Top.X, e.Bot.X);
	std::swap(e.Top.Z, e.Bot.Z);
}
//------------------------------------------------------------------------------

void SwapPoints(IntPoint& pt1, IntPoint& pt2);
//------------------------------------------------------------------------------

bool GetOverlapSegment(IntPoint pt1a, IntPoint pt1b, IntPoint pt2a,
	IntPoint pt2b, IntPoint& pt1, IntPoint& pt2);
//------------------------------------------------------------------------------

bool FirstIsBottomPt(const OutPt* btmPt1, const OutPt* btmPt2);
//------------------------------------------------------------------------------

OutPt* GetBottomPt(OutPt* pp);
//------------------------------------------------------------------------------
bool Pt2IsBetweenPt1AndPt3(const IntPoint pt1,
	const IntPoint pt2, const IntPoint pt3);
//------------------------------------------------------------------------------
bool HorzSegmentsOverlap(cInt seg1a, cInt seg1b, cInt seg2a, cInt seg2b);

DoublePoint GetUnitNormal(const IntPoint& pt1, const IntPoint& pt2);

//forward declarations (for stuff used internally) ...
struct TEdge;
struct IntersectNode;
struct LocalMinimum;
struct OutPt;
struct OutRec;
struct Join;

typedef std::vector < OutRec* > PolyOutList;
typedef std::vector < TEdge* > EdgeList;
typedef std::vector < Join* > JoinList;
typedef std::vector < IntersectNode* > IntersectList;

//------------------------------------------------------------------------------

//ClipperBase is the ancestor to the Clipper class. It should not be
//instantiated directly. This class simply abstracts the conversion of sets of
//polygon coordinates into edge objects that are stored in a LocalMinima list.
class ClipperBase
{
public:
  ClipperBase();
  virtual ~ClipperBase();
  virtual bool AddPath(const Path &pg, PolyType PolyTyp, bool Closed);
  bool AddPaths(const Paths &ppg, PolyType PolyTyp, bool Closed);
  virtual void Clear();
  IntRect GetBounds();
  bool PreserveCollinear() {return m_PreserveCollinear;};
  void PreserveCollinear(bool value) {m_PreserveCollinear = value;};
protected:
  void DisposeLocalMinimaList();
  TEdge* AddBoundsToLML(TEdge *e, bool IsClosed);
  virtual void Reset();
  TEdge* ProcessBound(TEdge* E, bool IsClockwise);
  void InsertScanbeam(const cInt Y);
  bool PopScanbeam(cInt &Y);
  bool LocalMinimaPending();
  bool PopLocalMinima(cInt Y, const LocalMinimum *&locMin);
  OutRec* CreateOutRec();
  void DisposeAllOutRecs();
  void DisposeOutRec(PolyOutList::size_type index);
  void SwapPositionsInAEL(TEdge *edge1, TEdge *edge2);
  void DeleteFromAEL(TEdge *e);
  void UpdateEdgeIntoAEL(TEdge *&e);

  typedef std::vector<LocalMinimum> MinimaList;
  MinimaList::iterator m_CurrentLM;
  MinimaList           m_MinimaList;

  bool              m_UseFullRange;
  EdgeList          m_edges;
  bool              m_PreserveCollinear;
  bool              m_HasOpenPaths;
  PolyOutList       m_PolyOuts;
  TEdge           *m_ActiveEdges;

  typedef std::priority_queue<cInt> ScanbeamList;
  ScanbeamList     m_Scanbeam;
};
//------------------------------------------------------------------------------

class Clipper : public virtual ClipperBase
{
public:
  Clipper(int initOptions = 0);
  bool Execute(ClipType clipType,
      Paths &solution,
      PolyFillType fillType = pftEvenOdd);
  bool Execute(ClipType clipType,
      Paths &solution,
      PolyFillType subjFillType,
      PolyFillType clipFillType);
  bool Execute(ClipType clipType,
      PolyTree &polytree,
      PolyFillType fillType = pftEvenOdd);
  bool Execute(ClipType clipType,
      PolyTree &polytree,
      PolyFillType subjFillType,
      PolyFillType clipFillType);
  bool ReverseSolution() { return m_ReverseOutput; };
  void ReverseSolution(bool value) {m_ReverseOutput = value;};
  bool StrictlySimple() {return m_StrictSimple;};
  void StrictlySimple(bool value) {m_StrictSimple = value;};
  //set the callback function for z value filling on intersections (otherwise Z is 0)
  void ZFillFunction(ZFillCallback zFillFunc);
protected:
  virtual bool ExecuteInternal();
private:
  JoinList         m_Joins;
  JoinList         m_GhostJoins;
  IntersectList    m_IntersectList;
  ClipType         m_ClipType;
  typedef std::list<cInt> MaximaList;
  MaximaList       m_Maxima;
  TEdge           *m_SortedEdges;
  bool             m_ExecuteLocked;
  PolyFillType     m_ClipFillType;
  PolyFillType     m_SubjFillType;
  bool             m_ReverseOutput;
  bool             m_UsingPolyTree; 
  bool             m_StrictSimple;
  ZFillCallback   m_ZFill; //custom callback 

  void SetWindingCount(TEdge& edge);
  bool IsEvenOddFillType(const TEdge& edge) const;
  bool IsEvenOddAltFillType(const TEdge& edge) const;
  void InsertLocalMinimaIntoAEL(const cInt botY);
  void InsertEdgeIntoAEL(TEdge *edge, TEdge* startEdge);
  void AddEdgeToSEL(TEdge *edge);
  bool PopEdgeFromSEL(TEdge *&edge);
  void CopyAELToSEL();
  void DeleteFromSEL(TEdge *e);
  void SwapPositionsInSEL(TEdge *edge1, TEdge *edge2);
  bool IsContributing(const TEdge& edge) const;
  bool IsTopHorz(const cInt XPos);
  void DoMaxima(TEdge *e);
  void ProcessHorizontals();
  void ProcessHorizontal(TEdge *horzEdge);
  void AddLocalMaxPoly(TEdge *e1, TEdge *e2, const IntPoint &pt);
  OutPt* AddLocalMinPoly(TEdge *e1, TEdge *e2, const IntPoint &pt);
  OutRec* GetOutRec(int idx);
  void AppendPolygon(TEdge *e1, TEdge *e2);
  void IntersectEdges(TEdge *e1, TEdge *e2, IntPoint &pt);
  OutPt* AddOutPt(TEdge *e, const IntPoint &pt);
  OutPt* GetLastOutPt(TEdge *e);
  bool ProcessIntersections(const cInt topY);
  void BuildIntersectList(const cInt topY);
  void ProcessIntersectList();
  void ProcessEdgesAtTopOfScanbeam(const cInt topY);
  void BuildResult(Paths& polys);
  void BuildResult2(PolyTree& polytree);
  void SetHoleState(TEdge *e, OutRec *outrec);
  void DisposeIntersectNodes();
  bool FixupIntersectionOrder();
  void FixupOutPolygon(OutRec &outrec);
  void FixupOutPolyline(OutRec &outrec);
  bool IsHole(TEdge *e);
  bool FindOwnerFromSplitRecs(OutRec &outRec, OutRec *&currOrfl);
  void FixHoleLinkage(OutRec &outrec);
  void AddJoin(OutPt *op1, OutPt *op2, const IntPoint offPt);
  void ClearJoins();
  void ClearGhostJoins();
  void AddGhostJoin(OutPt *op, const IntPoint offPt);
  bool JoinPoints(Join *j, OutRec* outRec1, OutRec* outRec2);
  void JoinCommonEdges();
  void DoSimplePolygons();
  void FixupFirstLefts1(OutRec* OldOutRec, OutRec* NewOutRec);
  void FixupFirstLefts2(OutRec* InnerOutRec, OutRec* OuterOutRec);
  void FixupFirstLefts3(OutRec* OldOutRec, OutRec* NewOutRec);
  void SetZ(IntPoint& pt, TEdge& e1, TEdge& e2);
};
//------------------------------------------------------------------------------

class ClipperOffset 
{
public:
  ClipperOffset(double miterLimit = 2.0, double roundPrecision = 0.25);
  ~ClipperOffset();
  void AddPath(const Path& path, JoinType joinType, EndType endType);
  void AddPaths(const Paths& paths, JoinType joinType, EndType endType);
  void Execute(Paths& solution, double delta);
  void Execute(PolyTree& solution, double delta);
  void ExecuteConst(Paths& solution, double delta, int step = -1);
  void Clear();
  double MiterLimit;
  double ArcTolerance;
private:
  Paths m_destPolys;
  Path m_srcPoly;
  Path m_destPoly;
  std::vector<DoublePoint> m_normals;
  double m_delta, m_sinA, m_sin, m_cos;
  double m_miterLim, m_StepsPerRad;
  IntPoint m_lowest;
  PolyNode m_polyNodes;

  void FixOrientations();
  void DoOffset(double delta, int step);
  void DoConstOffset(double delta, int step);
  void OffsetPoint(int j, int& k, JoinType jointype);
  void DoSquare(int j, int k);
  void DoMiter(int j, int k, double r);
  void DoRound(int j, int k);
};
//------------------------------------------------------------------------------

class clipperException : public std::exception
{
  public:
    clipperException(const char* description): m_descr(description) {}
    virtual ~clipperException() throw() {}
    virtual const char* what() const throw() {return m_descr.c_str();}
  private:
    std::string m_descr;
};
//------------------------------------------------------------------------------

class ClipperEx : public virtual ClipperBase
{
public:
	ClipperEx(int initOptions = 0);
	bool Execute(ClipType clipType,
		Paths& solution,
		PolyFillType fillType = pftEvenOdd);
	bool Execute(ClipType clipType,
		Paths& solution,
		PolyFillType subjFillType,
		PolyFillType clipFillType);
	bool Execute(ClipType clipType,
		PolyTree& polytree,
		PolyFillType fillType = pftEvenOdd);
	bool Execute(ClipType clipType,
		PolyTree& polytree,
		PolyFillType subjFillType,
		PolyFillType clipFillType);
	bool ReverseSolution() { return m_ReverseOutput; };
	void ReverseSolution(bool value) { m_ReverseOutput = value; };
	bool StrictlySimple() { return m_StrictSimple; };
	void StrictlySimple(bool value) { m_StrictSimple = value; };
	//set the callback function for z value filling on intersections (otherwise Z is 0)
	void ZFillFunction(ZFillCallback zFillFunc);
protected:
	virtual bool ExecuteInternal();
private:
	JoinList         m_Joins;
	JoinList         m_GhostJoins;
	IntersectList    m_IntersectList;
	ClipType         m_ClipType;
	typedef std::list<cInt> MaximaList;
	MaximaList       m_Maxima;
	TEdge* m_SortedEdges;
	bool             m_ExecuteLocked;
	PolyFillType     m_ClipFillType;
	PolyFillType     m_SubjFillType;
	bool             m_ReverseOutput;
	bool             m_UsingPolyTree;
	bool             m_StrictSimple;
	ZFillCallback   m_ZFill; //custom callback 
	void SetWindingCount(TEdge& edge);
	bool IsEvenOddFillType(const TEdge& edge) const;
	bool IsEvenOddAltFillType(const TEdge& edge) const;
	void InsertLocalMinimaIntoAEL(const cInt botY);
	void InsertEdgeIntoAEL(TEdge* edge, TEdge* startEdge);
	void AddEdgeToSEL(TEdge* edge);
	bool PopEdgeFromSEL(TEdge*& edge);
	void CopyAELToSEL();
	void DeleteFromSEL(TEdge* e);
	void SwapPositionsInSEL(TEdge* edge1, TEdge* edge2);
	bool IsContributing(const TEdge& edge) const;
	bool IsTopHorz(const cInt XPos);
	void DoMaxima(TEdge* e);
	void ProcessHorizontals();
	void ProcessHorizontal(TEdge* horzEdge);
	void AddLocalMaxPoly(TEdge* e1, TEdge* e2, const IntPoint& pt);
	OutPt* AddLocalMinPoly(TEdge* e1, TEdge* e2, const IntPoint& pt);
	OutRec* GetOutRec(int idx);
	void AppendPolygon(TEdge* e1, TEdge* e2);
	void IntersectEdges(TEdge* e1, TEdge* e2, IntPoint& pt);
	OutPt* AddOutPt(TEdge* e, const IntPoint& pt);
	OutPt* GetLastOutPt(TEdge* e);
	bool ProcessIntersections(const cInt topY);
	void BuildIntersectList(const cInt topY);
	void ProcessIntersectList();
	void ProcessEdgesAtTopOfScanbeam(const cInt topY);
	void BuildResult(Paths& polys);
	void BuildResult2(PolyTree& polytree);
	void SetHoleState(TEdge* e, OutRec* outrec);
	void DisposeIntersectNodes();
	bool FixupIntersectionOrder();
	void FixupOutPolygon(OutRec& outrec);
	void FixupOutPolyline(OutRec& outrec);
	bool IsHole(TEdge* e);
	bool FindOwnerFromSplitRecs(OutRec& outRec, OutRec*& currOrfl);
	void FixHoleLinkage(OutRec& outrec);
	void AddJoin(OutPt* op1, OutPt* op2, const IntPoint offPt);
	void ClearJoins();
	void ClearGhostJoins();
	void AddGhostJoin(OutPt* op, const IntPoint offPt);
	bool JoinPoints(Join* j, OutRec* outRec1, OutRec* outRec2);
	void JoinCommonEdges();
	void DoSimplePolygons();
	void FixupFirstLefts1(OutRec* OldOutRec, OutRec* NewOutRec);
	void FixupFirstLefts2(OutRec* InnerOutRec, OutRec* OuterOutRec);
	void FixupFirstLefts3(OutRec* OldOutRec, OutRec* NewOutRec);
	void SetZ(IntPoint& pt, TEdge& e1, TEdge& e2);
};

class ClipperOffsetEx
{
public:
	ClipperOffsetEx(double miterLimit = 2.0, double roundPrecision = 0.25);
	~ClipperOffsetEx();
	void AddPath(const Path& path, JoinType joinType, EndType endType);
	void AddPaths(const Paths& paths, JoinType joinType, EndType endType);
	void Execute(Paths& solution, double delta);
	void Execute(PolyTree& solution, double delta);
	void ExecuteConst(Paths& solution, double delta, int step = -1);
	void Clear();
	double MiterLimit;
	double ArcTolerance;
private:
	Paths m_destPolys;
	Path m_srcPoly;
	Path m_destPoly;
	std::vector<DoublePoint> m_normals;
	double m_delta, m_sinA, m_sin, m_cos;
	double m_miterLim, m_StepsPerRad;
	IntPoint m_lowest;
	PolyNode m_polyNodes;

	void FixOrientations();
	void DoOffset(double delta, int step);
	void DoConstOffset(double delta, int step);
	void OffsetPoint(int j, int& k, JoinType jointype);
	void DoSquare(int j, int k);
	void DoMiter(int j, int k, double r);
	void DoRound(int j, int k);
};

} //ClipperLibXYZ namespace

namespace std {
	template <>
	struct hash<ClipperLibXYZ::IntPoint> {
		size_t operator()(const ClipperLibXYZ::IntPoint& pp) const
		{
			static int prime = 31;
			int result = 89;
			result = result * prime + pp.X;
			result = result * prime + pp.Y;
			return result;
		}
	};
}

#define INT2MM(n) (float(n) / 1000.0f)
typedef std::function<void(ClipperLibXYZ::PolyNode*)> polyNodeFunc;
typedef std::function<void(std::vector<ClipperLibXYZ::PolyNode*>&)> polyLevelFunc;

#endif //clipper_hpp_xyz


