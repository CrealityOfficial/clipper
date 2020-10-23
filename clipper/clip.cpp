#include <clipper.hpp>

namespace ClipperLib
{
	ClipperEx::ClipperEx(int initOptions) : ClipperBase() //constructor
	{
		m_ExecuteLocked = false;
		m_UseFullRange = false;
		m_ReverseOutput = ((initOptions & ioReverseSolution) != 0);
		m_StrictSimple = ((initOptions & ioStrictlySimple) != 0);
		m_PreserveCollinear = ((initOptions & ioPreserveCollinear) != 0);
		m_HasOpenPaths = false;
#ifdef use_xyz  
		m_ZFill = 0;
#endif
	}
	//------------------------------------------------------------------------------

#ifdef use_xyz  
	void ClipperEx::ZFillFunction(ZFillCallback zFillFunc)
	{
		m_ZFill = zFillFunc;
	}
	//------------------------------------------------------------------------------
#endif

	bool ClipperEx::Execute(ClipType clipType, Paths& solution, PolyFillType fillType)
	{
		return Execute(clipType, solution, fillType, fillType);
	}
	//------------------------------------------------------------------------------

	bool ClipperEx::Execute(ClipType clipType, PolyTree& polytree, PolyFillType fillType)
	{
		return Execute(clipType, polytree, fillType, fillType);
	}
	//------------------------------------------------------------------------------

	bool ClipperEx::Execute(ClipType clipType, Paths& solution,
		PolyFillType subjFillType, PolyFillType clipFillType)
	{
		if (m_ExecuteLocked) return false;
		if (m_HasOpenPaths)
			throw clipperException("Error: PolyTree struct is needed for open path clipping.");
		m_ExecuteLocked = true;
		solution.resize(0);
		m_SubjFillType = subjFillType;
		m_ClipFillType = clipFillType;
		m_ClipType = clipType;
		m_UsingPolyTree = false;
		bool succeeded = ExecuteInternal();
		if (succeeded) BuildResult(solution);
		DisposeAllOutRecs();
		m_ExecuteLocked = false;
		return succeeded;
	}
	//------------------------------------------------------------------------------

	bool ClipperEx::Execute(ClipType clipType, PolyTree& polytree,
		PolyFillType subjFillType, PolyFillType clipFillType)
	{
		if (m_ExecuteLocked) return false;
		m_ExecuteLocked = true;
		m_SubjFillType = subjFillType;
		m_ClipFillType = clipFillType;
		m_ClipType = clipType;
		m_UsingPolyTree = true;
		bool succeeded = ExecuteInternal();
		if (succeeded) BuildResult2(polytree);
		DisposeAllOutRecs();
		m_ExecuteLocked = false;
		return succeeded;
	}
	//------------------------------------------------------------------------------

	void ClipperEx::FixHoleLinkage(OutRec& outrec)
	{
		//skip OutRecs that (a) contain outermost polygons or
		//(b) already have the correct owner/child linkage ...
		if (!outrec.FirstLeft ||
			(outrec.IsHole != outrec.FirstLeft->IsHole &&
				outrec.FirstLeft->Pts)) return;

		OutRec* orfl = outrec.FirstLeft;
		while (orfl && ((orfl->IsHole == outrec.IsHole) || !orfl->Pts))
			orfl = orfl->FirstLeft;
		outrec.FirstLeft = orfl;
	}
	//------------------------------------------------------------------------------

	bool ClipperEx::ExecuteInternal()
	{
		bool succeeded = true;
		try {
			Reset();
			m_Maxima = MaximaList();
			m_SortedEdges = 0;

			succeeded = true;
			cInt botY, topY;
			if (!PopScanbeam(botY)) return false;
			InsertLocalMinimaIntoAEL(botY);
			while (PopScanbeam(topY) || LocalMinimaPending())
			{
				ProcessHorizontals();
				ClearGhostJoins();
				if (!ProcessIntersections(topY))
				{
					succeeded = false;
					break;
				}
				ProcessEdgesAtTopOfScanbeam(topY);
				botY = topY;
				InsertLocalMinimaIntoAEL(botY);
			}
		}
		catch (...)
		{
			succeeded = false;
		}

		if (succeeded)
		{
			//fix orientations ...
			for (PolyOutList::size_type i = 0; i < m_PolyOuts.size(); ++i)
			{
				OutRec* outRec = m_PolyOuts[i];
				if (!outRec->Pts || outRec->IsOpen) continue;
				if ((outRec->IsHole ^ m_ReverseOutput) == (Area(*outRec) > 0))
					ReversePolyPtLinks(outRec->Pts);
			}

			if (!m_Joins.empty()) JoinCommonEdges();

			//unfortunately FixupOutPolygon() must be done after JoinCommonEdges()
			for (PolyOutList::size_type i = 0; i < m_PolyOuts.size(); ++i)
			{
				OutRec* outRec = m_PolyOuts[i];
				if (!outRec->Pts) continue;
				if (outRec->IsOpen)
					FixupOutPolyline(*outRec);
				else
					FixupOutPolygon(*outRec);
			}

			if (m_StrictSimple) DoSimplePolygons();
		}

		ClearJoins();
		ClearGhostJoins();
		return succeeded;
	}
	//------------------------------------------------------------------------------

	void ClipperEx::SetWindingCount(TEdge& edge)
	{
		TEdge* e = edge.PrevInAEL;
		//find the edge of the same polytype that immediately preceeds 'edge' in AEL
		while (e && ((e->PolyTyp != edge.PolyTyp) || (e->WindDelta == 0))) e = e->PrevInAEL;
		if (!e)
		{
			if (edge.WindDelta == 0)
			{
				PolyFillType pft = (edge.PolyTyp == ptSubject ? m_SubjFillType : m_ClipFillType);
				edge.WindCnt = (pft == pftNegative ? -1 : 1);
			}
			else
				edge.WindCnt = edge.WindDelta;
			edge.WindCnt2 = 0;
			e = m_ActiveEdges; //ie get ready to calc WindCnt2
		}
		else if (edge.WindDelta == 0 && m_ClipType != ctUnion)
		{
			edge.WindCnt = 1;
			edge.WindCnt2 = e->WindCnt2;
			e = e->NextInAEL; //ie get ready to calc WindCnt2
		}
		else if (IsEvenOddFillType(edge))
		{
			//EvenOdd filling ...
			if (edge.WindDelta == 0)
			{
				//are we inside a subj polygon ...
				bool Inside = true;
				TEdge* e2 = e->PrevInAEL;
				while (e2)
				{
					if (e2->PolyTyp == e->PolyTyp && e2->WindDelta != 0)
						Inside = !Inside;
					e2 = e2->PrevInAEL;
				}
				edge.WindCnt = (Inside ? 0 : 1);
			}
			else
			{
				edge.WindCnt = edge.WindDelta;
			}
			edge.WindCnt2 = e->WindCnt2;
			e = e->NextInAEL; //ie get ready to calc WindCnt2
		}
		else
		{
			//nonZero, Positive or Negative filling ...
			if (e->WindCnt * e->WindDelta < 0)
			{
				//prev edge is 'decreasing' WindCount (WC) toward zero
				//so we're outside the previous polygon ...
				if (Abs(e->WindCnt) > 1)
				{
					//outside prev poly but still inside another.
					//when reversing direction of prev poly use the same WC 
					if (e->WindDelta * edge.WindDelta < 0) edge.WindCnt = e->WindCnt;
					//otherwise continue to 'decrease' WC ...
					else edge.WindCnt = e->WindCnt + edge.WindDelta;
				}
				else
					//now outside all polys of same polytype so set own WC ...
					edge.WindCnt = (edge.WindDelta == 0 ? 1 : edge.WindDelta);
			}
			else
			{
				//prev edge is 'increasing' WindCount (WC) away from zero
				//so we're inside the previous polygon ...
				if (edge.WindDelta == 0)
					edge.WindCnt = (e->WindCnt < 0 ? e->WindCnt - 1 : e->WindCnt + 1);
				//if wind direction is reversing prev then use same WC
				else if (e->WindDelta * edge.WindDelta < 0) edge.WindCnt = e->WindCnt;
				//otherwise add to WC ...
				else edge.WindCnt = e->WindCnt + edge.WindDelta;
			}
			edge.WindCnt2 = e->WindCnt2;
			e = e->NextInAEL; //ie get ready to calc WindCnt2
		}

		//update WindCnt2 ...
		if (IsEvenOddAltFillType(edge))
		{
			//EvenOdd filling ...
			while (e != &edge)
			{
				if (e->WindDelta != 0)
					edge.WindCnt2 = (edge.WindCnt2 == 0 ? 1 : 0);
				e = e->NextInAEL;
			}
		}
		else
		{
			//nonZero, Positive or Negative filling ...
			while (e != &edge)
			{
				edge.WindCnt2 += e->WindDelta;
				e = e->NextInAEL;
			}
		}
	}
	//------------------------------------------------------------------------------

	bool ClipperEx::IsEvenOddFillType(const TEdge& edge) const
	{
		if (edge.PolyTyp == ptSubject)
			return m_SubjFillType == pftEvenOdd; else
			return m_ClipFillType == pftEvenOdd;
	}
	//------------------------------------------------------------------------------

	bool ClipperEx::IsEvenOddAltFillType(const TEdge& edge) const
	{
		if (edge.PolyTyp == ptSubject)
			return m_ClipFillType == pftEvenOdd; else
			return m_SubjFillType == pftEvenOdd;
	}
	//------------------------------------------------------------------------------

	bool ClipperEx::IsContributing(const TEdge& edge) const
	{
		PolyFillType pft, pft2;
		if (edge.PolyTyp == ptSubject)
		{
			pft = m_SubjFillType;
			pft2 = m_ClipFillType;
		}
		else
		{
			pft = m_ClipFillType;
			pft2 = m_SubjFillType;
		}

		switch (pft)
		{
		case pftEvenOdd:
			//return false if a subj line has been flagged as inside a subj polygon
			if (edge.WindDelta == 0 && edge.WindCnt != 1) return false;
			break;
		case pftNonZero:
			if (Abs(edge.WindCnt) != 1) return false;
			break;
		case pftPositive:
			if (edge.WindCnt != 1) return false;
			break;
		default: //pftNegative
			if (edge.WindCnt != -1) return false;
		}

		switch (m_ClipType)
		{
		case ctIntersection:
			switch (pft2)
			{
			case pftEvenOdd:
			case pftNonZero:
				return (edge.WindCnt2 != 0);
			case pftPositive:
				return (edge.WindCnt2 > 0);
			default:
				return (edge.WindCnt2 < 0);
			}
			break;
		case ctUnion:
			switch (pft2)
			{
			case pftEvenOdd:
			case pftNonZero:
				return (edge.WindCnt2 == 0);
			case pftPositive:
				return (edge.WindCnt2 <= 0);
			default:
				return (edge.WindCnt2 >= 0);
			}
			break;
		case ctDifference:
			if (edge.PolyTyp == ptSubject)
				switch (pft2)
				{
				case pftEvenOdd:
				case pftNonZero:
					return (edge.WindCnt2 == 0);
				case pftPositive:
					return (edge.WindCnt2 <= 0);
				default:
					return (edge.WindCnt2 >= 0);
				}
			else
				switch (pft2)
				{
				case pftEvenOdd:
				case pftNonZero:
					return (edge.WindCnt2 != 0);
				case pftPositive:
					return (edge.WindCnt2 > 0);
				default:
					return (edge.WindCnt2 < 0);
				}
			break;
		case ctXor:
			if (edge.WindDelta == 0) //XOr always contributing unless open
				switch (pft2)
				{
				case pftEvenOdd:
				case pftNonZero:
					return (edge.WindCnt2 == 0);
				case pftPositive:
					return (edge.WindCnt2 <= 0);
				default:
					return (edge.WindCnt2 >= 0);
				}
			else
				return true;
			break;
		default:
			return true;
		}
	}
	//------------------------------------------------------------------------------

	OutPt* ClipperEx::AddLocalMinPoly(TEdge* e1, TEdge* e2, const IntPoint& Pt)
	{
		OutPt* result;
		TEdge* e, * prevE;
		if (IsHorizontal(*e2) || (e1->Dx > e2->Dx))
		{
			result = AddOutPt(e1, Pt);
			e2->OutIdx = e1->OutIdx;
			e1->Side = esLeft;
			e2->Side = esRight;
			e = e1;
			if (e->PrevInAEL == e2)
				prevE = e2->PrevInAEL;
			else
				prevE = e->PrevInAEL;
		}
		else
		{
			result = AddOutPt(e2, Pt);
			e1->OutIdx = e2->OutIdx;
			e1->Side = esRight;
			e2->Side = esLeft;
			e = e2;
			if (e->PrevInAEL == e1)
				prevE = e1->PrevInAEL;
			else
				prevE = e->PrevInAEL;
		}

		if (prevE && prevE->OutIdx >= 0 && prevE->Top.Y < Pt.Y && e->Top.Y < Pt.Y)
		{
			cInt xPrev = TopX(*prevE, Pt.Y);
			cInt xE = TopX(*e, Pt.Y);
			if (xPrev == xE && (e->WindDelta != 0) && (prevE->WindDelta != 0) &&
				SlopesEqual(IntPoint(xPrev, Pt.Y), prevE->Top, IntPoint(xE, Pt.Y), e->Top, m_UseFullRange))
			{
				OutPt* outPt = AddOutPt(prevE, Pt);
				AddJoin(result, outPt, e->Top);
			}
		}
		return result;
	}
	//------------------------------------------------------------------------------

	void ClipperEx::AddLocalMaxPoly(TEdge* e1, TEdge* e2, const IntPoint& Pt)
	{
		AddOutPt(e1, Pt);
		if (e2->WindDelta == 0) AddOutPt(e2, Pt);
		if (e1->OutIdx == e2->OutIdx)
		{
			e1->OutIdx = Unassigned;
			e2->OutIdx = Unassigned;
		}
		else if (e1->OutIdx < e2->OutIdx)
			AppendPolygon(e1, e2);
		else
			AppendPolygon(e2, e1);
	}
	//------------------------------------------------------------------------------

	void ClipperEx::AddEdgeToSEL(TEdge* edge)
	{
		//SEL pointers in PEdge are reused to build a list of horizontal edges.
		//However, we don't need to worry about order with horizontal edge processing.
		if (!m_SortedEdges)
		{
			m_SortedEdges = edge;
			edge->PrevInSEL = 0;
			edge->NextInSEL = 0;
		}
		else
		{
			edge->NextInSEL = m_SortedEdges;
			edge->PrevInSEL = 0;
			m_SortedEdges->PrevInSEL = edge;
			m_SortedEdges = edge;
		}
	}
	//------------------------------------------------------------------------------

	bool ClipperEx::PopEdgeFromSEL(TEdge*& edge)
	{
		if (!m_SortedEdges) return false;
		edge = m_SortedEdges;
		DeleteFromSEL(m_SortedEdges);
		return true;
	}
	//------------------------------------------------------------------------------

	void ClipperEx::CopyAELToSEL()
	{
		TEdge* e = m_ActiveEdges;
		m_SortedEdges = e;
		while (e)
		{
			e->PrevInSEL = e->PrevInAEL;
			e->NextInSEL = e->NextInAEL;
			e = e->NextInAEL;
		}
	}
	//------------------------------------------------------------------------------

	void ClipperEx::AddJoin(OutPt* op1, OutPt* op2, const IntPoint OffPt)
	{
		Join* j = new Join;
		j->OutPt1 = op1;
		j->OutPt2 = op2;
		j->OffPt = OffPt;
		m_Joins.push_back(j);
	}
	//------------------------------------------------------------------------------

	void ClipperEx::ClearJoins()
	{
		for (JoinList::size_type i = 0; i < m_Joins.size(); i++)
			delete m_Joins[i];
		m_Joins.resize(0);
	}
	//------------------------------------------------------------------------------

	void ClipperEx::ClearGhostJoins()
	{
		for (JoinList::size_type i = 0; i < m_GhostJoins.size(); i++)
			delete m_GhostJoins[i];
		m_GhostJoins.resize(0);
	}
	//------------------------------------------------------------------------------

	void ClipperEx::AddGhostJoin(OutPt* op, const IntPoint OffPt)
	{
		Join* j = new Join;
		j->OutPt1 = op;
		j->OutPt2 = 0;
		j->OffPt = OffPt;
		m_GhostJoins.push_back(j);
	}
	//------------------------------------------------------------------------------

	void ClipperEx::InsertLocalMinimaIntoAEL(const cInt botY)
	{
		const LocalMinimum* lm;
		while (PopLocalMinima(botY, lm))
		{
			TEdge* lb = lm->LeftBound;
			TEdge* rb = lm->RightBound;

			OutPt* Op1 = 0;
			if (!lb)
			{
				//nb: don't insert LB into either AEL or SEL
				InsertEdgeIntoAEL(rb, 0);
				SetWindingCount(*rb);
				if (IsContributing(*rb))
					Op1 = AddOutPt(rb, rb->Bot);
			}
			else if (!rb)
			{
				InsertEdgeIntoAEL(lb, 0);
				SetWindingCount(*lb);
				if (IsContributing(*lb))
					Op1 = AddOutPt(lb, lb->Bot);
				InsertScanbeam(lb->Top.Y);
			}
			else
			{
				InsertEdgeIntoAEL(lb, 0);
				InsertEdgeIntoAEL(rb, lb);
				SetWindingCount(*lb);
				rb->WindCnt = lb->WindCnt;
				rb->WindCnt2 = lb->WindCnt2;
				if (IsContributing(*lb))
					Op1 = AddLocalMinPoly(lb, rb, lb->Bot);
				InsertScanbeam(lb->Top.Y);
			}

			if (rb)
			{
				if (IsHorizontal(*rb))
				{
					AddEdgeToSEL(rb);
					if (rb->NextInLML)
						InsertScanbeam(rb->NextInLML->Top.Y);
				}
				else InsertScanbeam(rb->Top.Y);
			}

			if (!lb || !rb) continue;

			//if any output polygons share an edge, they'll need joining later ...
			if (Op1 && IsHorizontal(*rb) &&
				m_GhostJoins.size() > 0 && (rb->WindDelta != 0))
			{
				for (JoinList::size_type i = 0; i < m_GhostJoins.size(); ++i)
				{
					Join* jr = m_GhostJoins[i];
					//if the horizontal Rb and a 'ghost' horizontal overlap, then convert
					//the 'ghost' join to a real join ready for later ...
					if (HorzSegmentsOverlap(jr->OutPt1->Pt.X, jr->OffPt.X, rb->Bot.X, rb->Top.X))
						AddJoin(jr->OutPt1, Op1, jr->OffPt);
				}
			}

			if (lb->OutIdx >= 0 && lb->PrevInAEL &&
				lb->PrevInAEL->Curr.X == lb->Bot.X &&
				lb->PrevInAEL->OutIdx >= 0 &&
				SlopesEqual(lb->PrevInAEL->Bot, lb->PrevInAEL->Top, lb->Curr, lb->Top, m_UseFullRange) &&
				(lb->WindDelta != 0) && (lb->PrevInAEL->WindDelta != 0))
			{
				OutPt* Op2 = AddOutPt(lb->PrevInAEL, lb->Bot);
				AddJoin(Op1, Op2, lb->Top);
			}

			if (lb->NextInAEL != rb)
			{

				if (rb->OutIdx >= 0 && rb->PrevInAEL->OutIdx >= 0 &&
					SlopesEqual(rb->PrevInAEL->Curr, rb->PrevInAEL->Top, rb->Curr, rb->Top, m_UseFullRange) &&
					(rb->WindDelta != 0) && (rb->PrevInAEL->WindDelta != 0))
				{
					OutPt* Op2 = AddOutPt(rb->PrevInAEL, rb->Bot);
					AddJoin(Op1, Op2, rb->Top);
				}

				TEdge* e = lb->NextInAEL;
				if (e)
				{
					while (e != rb)
					{
						//nb: For calculating winding counts etc, IntersectEdges() assumes
						//that param1 will be to the Right of param2 ABOVE the intersection ...
						IntersectEdges(rb, e, lb->Curr); //order important here
						e = e->NextInAEL;
					}
				}
			}

		}
	}
	//------------------------------------------------------------------------------

	void ClipperEx::DeleteFromSEL(TEdge* e)
	{
		TEdge* SelPrev = e->PrevInSEL;
		TEdge* SelNext = e->NextInSEL;
		if (!SelPrev && !SelNext && (e != m_SortedEdges)) return; //already deleted
		if (SelPrev) SelPrev->NextInSEL = SelNext;
		else m_SortedEdges = SelNext;
		if (SelNext) SelNext->PrevInSEL = SelPrev;
		e->NextInSEL = 0;
		e->PrevInSEL = 0;
	}
	//------------------------------------------------------------------------------

#ifdef use_xyz
	void ClipperEx::SetZ(IntPoint& pt, TEdge& e1, TEdge& e2)
	{
		if (pt.Z != 0 || !m_ZFill) return;
		else if (pt == e1.Bot) pt.Z = e1.Bot.Z;
		else if (pt == e1.Top) pt.Z = e1.Top.Z;
		else if (pt == e2.Bot) pt.Z = e2.Bot.Z;
		else if (pt == e2.Top) pt.Z = e2.Top.Z;
		else (*m_ZFill)(e1.Bot, e1.Top, e2.Bot, e2.Top, pt);
	}
	//------------------------------------------------------------------------------
#endif

	void ClipperEx::IntersectEdges(TEdge* e1, TEdge* e2, IntPoint& Pt)
	{
		bool e1Contributing = (e1->OutIdx >= 0);
		bool e2Contributing = (e2->OutIdx >= 0);

#ifdef use_xyz
		SetZ(Pt, *e1, *e2);
#endif

#ifdef use_lines
		//if either edge is on an OPEN path ...
		if (e1->WindDelta == 0 || e2->WindDelta == 0)
		{
			//ignore subject-subject open path intersections UNLESS they
			//are both open paths, AND they are both 'contributing maximas' ...
			if (e1->WindDelta == 0 && e2->WindDelta == 0) return;

			//if intersecting a subj line with a subj poly ...
			else if (e1->PolyTyp == e2->PolyTyp &&
				e1->WindDelta != e2->WindDelta && m_ClipType == ctUnion)
			{
				if (e1->WindDelta == 0)
				{
					if (e2Contributing)
					{
						AddOutPt(e1, Pt);
						if (e1Contributing) e1->OutIdx = Unassigned;
					}
				}
				else
				{
					if (e1Contributing)
					{
						AddOutPt(e2, Pt);
						if (e2Contributing) e2->OutIdx = Unassigned;
					}
				}
			}
			else if (e1->PolyTyp != e2->PolyTyp)
			{
				//toggle subj open path OutIdx on/off when Abs(clip.WndCnt) == 1 ...
				if ((e1->WindDelta == 0) && abs(e2->WindCnt) == 1 &&
					(m_ClipType != ctUnion || e2->WindCnt2 == 0))
				{
					AddOutPt(e1, Pt);
					if (e1Contributing) e1->OutIdx = Unassigned;
				}
				else if ((e2->WindDelta == 0) && (abs(e1->WindCnt) == 1) &&
					(m_ClipType != ctUnion || e1->WindCnt2 == 0))
				{
					AddOutPt(e2, Pt);
					if (e2Contributing) e2->OutIdx = Unassigned;
				}
			}
			return;
		}
#endif

		//update winding counts...
		//assumes that e1 will be to the Right of e2 ABOVE the intersection
		if (e1->PolyTyp == e2->PolyTyp)
		{
			if (IsEvenOddFillType(*e1))
			{
				int oldE1WindCnt = e1->WindCnt;
				e1->WindCnt = e2->WindCnt;
				e2->WindCnt = oldE1WindCnt;
			}
			else
			{
				if (e1->WindCnt + e2->WindDelta == 0) e1->WindCnt = -e1->WindCnt;
				else e1->WindCnt += e2->WindDelta;
				if (e2->WindCnt - e1->WindDelta == 0) e2->WindCnt = -e2->WindCnt;
				else e2->WindCnt -= e1->WindDelta;
			}
		}
		else
		{
			if (!IsEvenOddFillType(*e2)) e1->WindCnt2 += e2->WindDelta;
			else e1->WindCnt2 = (e1->WindCnt2 == 0) ? 1 : 0;
			if (!IsEvenOddFillType(*e1)) e2->WindCnt2 -= e1->WindDelta;
			else e2->WindCnt2 = (e2->WindCnt2 == 0) ? 1 : 0;
		}

		PolyFillType e1FillType, e2FillType, e1FillType2, e2FillType2;
		if (e1->PolyTyp == ptSubject)
		{
			e1FillType = m_SubjFillType;
			e1FillType2 = m_ClipFillType;
		}
		else
		{
			e1FillType = m_ClipFillType;
			e1FillType2 = m_SubjFillType;
		}
		if (e2->PolyTyp == ptSubject)
		{
			e2FillType = m_SubjFillType;
			e2FillType2 = m_ClipFillType;
		}
		else
		{
			e2FillType = m_ClipFillType;
			e2FillType2 = m_SubjFillType;
		}

		cInt e1Wc, e2Wc;
		switch (e1FillType)
		{
		case pftPositive: e1Wc = e1->WindCnt; break;
		case pftNegative: e1Wc = -e1->WindCnt; break;
		default: e1Wc = Abs(e1->WindCnt);
		}
		switch (e2FillType)
		{
		case pftPositive: e2Wc = e2->WindCnt; break;
		case pftNegative: e2Wc = -e2->WindCnt; break;
		default: e2Wc = Abs(e2->WindCnt);
		}

		if (e1Contributing && e2Contributing)
		{
			if ((e1Wc != 0 && e1Wc != 1) || (e2Wc != 0 && e2Wc != 1) ||
				(e1->PolyTyp != e2->PolyTyp && m_ClipType != ctXor))
			{
				AddLocalMaxPoly(e1, e2, Pt);
			}
			else
			{
				AddOutPt(e1, Pt);
				AddOutPt(e2, Pt);
				SwapSides(*e1, *e2);
				SwapPolyIndexes(*e1, *e2);
			}
		}
		else if (e1Contributing)
		{
			if (e2Wc == 0 || e2Wc == 1)
			{
				AddOutPt(e1, Pt);
				SwapSides(*e1, *e2);
				SwapPolyIndexes(*e1, *e2);
			}
		}
		else if (e2Contributing)
		{
			if (e1Wc == 0 || e1Wc == 1)
			{
				AddOutPt(e2, Pt);
				SwapSides(*e1, *e2);
				SwapPolyIndexes(*e1, *e2);
			}
		}
		else if ((e1Wc == 0 || e1Wc == 1) && (e2Wc == 0 || e2Wc == 1))
		{
			//neither edge is currently contributing ...

			cInt e1Wc2, e2Wc2;
			switch (e1FillType2)
			{
			case pftPositive: e1Wc2 = e1->WindCnt2; break;
			case pftNegative: e1Wc2 = -e1->WindCnt2; break;
			default: e1Wc2 = Abs(e1->WindCnt2);
			}
			switch (e2FillType2)
			{
			case pftPositive: e2Wc2 = e2->WindCnt2; break;
			case pftNegative: e2Wc2 = -e2->WindCnt2; break;
			default: e2Wc2 = Abs(e2->WindCnt2);
			}

			if (e1->PolyTyp != e2->PolyTyp)
			{
				AddLocalMinPoly(e1, e2, Pt);
			}
			else if (e1Wc == 1 && e2Wc == 1)
				switch (m_ClipType) {
				case ctIntersection:
					if (e1Wc2 > 0 && e2Wc2 > 0)
						AddLocalMinPoly(e1, e2, Pt);
					break;
				case ctUnion:
					if (e1Wc2 <= 0 && e2Wc2 <= 0)
						AddLocalMinPoly(e1, e2, Pt);
					break;
				case ctDifference:
					if (((e1->PolyTyp == ptClip) && (e1Wc2 > 0) && (e2Wc2 > 0)) ||
						((e1->PolyTyp == ptSubject) && (e1Wc2 <= 0) && (e2Wc2 <= 0)))
						AddLocalMinPoly(e1, e2, Pt);
					break;
				case ctXor:
					AddLocalMinPoly(e1, e2, Pt);
				}
			else
				SwapSides(*e1, *e2);
		}
	}
	//------------------------------------------------------------------------------

	void ClipperEx::SetHoleState(TEdge* e, OutRec* outrec)
	{
		TEdge* e2 = e->PrevInAEL;
		TEdge* eTmp = 0;
		while (e2)
		{
			if (e2->OutIdx >= 0 && e2->WindDelta != 0)
			{
				if (!eTmp) eTmp = e2;
				else if (eTmp->OutIdx == e2->OutIdx) eTmp = 0;
			}
			e2 = e2->PrevInAEL;
		}
		if (!eTmp)
		{
			outrec->FirstLeft = 0;
			outrec->IsHole = false;
		}
		else
		{
			outrec->FirstLeft = m_PolyOuts[eTmp->OutIdx];
			outrec->IsHole = !outrec->FirstLeft->IsHole;
		}
	}
	//------------------------------------------------------------------------------
	//------------------------------------------------------------------------------

	OutRec* ClipperEx::GetOutRec(int Idx)
	{
		OutRec* outrec = m_PolyOuts[Idx];
		while (outrec != m_PolyOuts[outrec->Idx])
			outrec = m_PolyOuts[outrec->Idx];
		return outrec;
	}
	//------------------------------------------------------------------------------

	void ClipperEx::AppendPolygon(TEdge* e1, TEdge* e2)
	{
		//get the start and ends of both output polygons ...
		OutRec* outRec1 = m_PolyOuts[e1->OutIdx];
		OutRec* outRec2 = m_PolyOuts[e2->OutIdx];

		OutRec* holeStateRec;
		if (OutRec1RightOfOutRec2(outRec1, outRec2))
			holeStateRec = outRec2;
		else if (OutRec1RightOfOutRec2(outRec2, outRec1))
			holeStateRec = outRec1;
		else
			holeStateRec = GetLowermostRec(outRec1, outRec2);

		//get the start and ends of both output polygons and
		//join e2 poly onto e1 poly and delete pointers to e2 ...

		OutPt* p1_lft = outRec1->Pts;
		OutPt* p1_rt = p1_lft->Prev;
		OutPt* p2_lft = outRec2->Pts;
		OutPt* p2_rt = p2_lft->Prev;

		//join e2 poly onto e1 poly and delete pointers to e2 ...
		if (e1->Side == esLeft)
		{
			if (e2->Side == esLeft)
			{
				//z y x a b c
				ReversePolyPtLinks(p2_lft);
				p2_lft->Next = p1_lft;
				p1_lft->Prev = p2_lft;
				p1_rt->Next = p2_rt;
				p2_rt->Prev = p1_rt;
				outRec1->Pts = p2_rt;
			}
			else
			{
				//x y z a b c
				p2_rt->Next = p1_lft;
				p1_lft->Prev = p2_rt;
				p2_lft->Prev = p1_rt;
				p1_rt->Next = p2_lft;
				outRec1->Pts = p2_lft;
			}
		}
		else
		{
			if (e2->Side == esRight)
			{
				//a b c z y x
				ReversePolyPtLinks(p2_lft);
				p1_rt->Next = p2_rt;
				p2_rt->Prev = p1_rt;
				p2_lft->Next = p1_lft;
				p1_lft->Prev = p2_lft;
			}
			else
			{
				//a b c x y z
				p1_rt->Next = p2_lft;
				p2_lft->Prev = p1_rt;
				p1_lft->Prev = p2_rt;
				p2_rt->Next = p1_lft;
			}
		}

		outRec1->BottomPt = 0;
		if (holeStateRec == outRec2)
		{
			if (outRec2->FirstLeft != outRec1)
				outRec1->FirstLeft = outRec2->FirstLeft;
			outRec1->IsHole = outRec2->IsHole;
		}
		outRec2->Pts = 0;
		outRec2->BottomPt = 0;
		outRec2->FirstLeft = outRec1;

		int OKIdx = e1->OutIdx;
		int ObsoleteIdx = e2->OutIdx;

		e1->OutIdx = Unassigned; //nb: safe because we only get here via AddLocalMaxPoly
		e2->OutIdx = Unassigned;

		TEdge* e = m_ActiveEdges;
		while (e)
		{
			if (e->OutIdx == ObsoleteIdx)
			{
				e->OutIdx = OKIdx;
				e->Side = e1->Side;
				break;
			}
			e = e->NextInAEL;
		}

		outRec2->Idx = outRec1->Idx;
	}
	//------------------------------------------------------------------------------

	OutPt* ClipperEx::AddOutPt(TEdge* e, const IntPoint& pt)
	{
		if (e->OutIdx < 0)
		{
			OutRec* outRec = CreateOutRec();
			outRec->IsOpen = (e->WindDelta == 0);
			OutPt* newOp = new OutPt;
			outRec->Pts = newOp;
			newOp->Idx = outRec->Idx;
			newOp->Pt = pt;
			newOp->Next = newOp;
			newOp->Prev = newOp;
			if (!outRec->IsOpen)
				SetHoleState(e, outRec);
			e->OutIdx = outRec->Idx;
			return newOp;
		}
		else
		{
			OutRec* outRec = m_PolyOuts[e->OutIdx];
			//OutRec.Pts is the 'Left-most' point & OutRec.Pts.Prev is the 'Right-most'
			OutPt* op = outRec->Pts;

			bool ToFront = (e->Side == esLeft);
			if (ToFront && (pt == op->Pt)) return op;
			else if (!ToFront && (pt == op->Prev->Pt)) return op->Prev;

			OutPt* newOp = new OutPt;
			newOp->Idx = outRec->Idx;
			newOp->Pt = pt;
			newOp->Next = op;
			newOp->Prev = op->Prev;
			newOp->Prev->Next = newOp;
			op->Prev = newOp;
			if (ToFront) outRec->Pts = newOp;
			return newOp;
		}
	}
	//------------------------------------------------------------------------------

	OutPt* ClipperEx::GetLastOutPt(TEdge* e)
	{
		OutRec* outRec = m_PolyOuts[e->OutIdx];
		if (e->Side == esLeft)
			return outRec->Pts;
		else
			return outRec->Pts->Prev;
	}
	//------------------------------------------------------------------------------

	void ClipperEx::ProcessHorizontals()
	{
		TEdge* horzEdge;
		while (PopEdgeFromSEL(horzEdge))
			ProcessHorizontal(horzEdge);
	}
	//------------------------------------------------------------------------------

	inline bool IsMinima(TEdge* e)
	{
		return e && (e->Prev->NextInLML != e) && (e->Next->NextInLML != e);
	}
	//------------------------------------------------------------------------------

	inline bool IsMaxima(TEdge* e, const cInt Y)
	{
		return e && e->Top.Y == Y && !e->NextInLML;
	}
	//------------------------------------------------------------------------------

	inline bool IsIntermediate(TEdge* e, const cInt Y)
	{
		return e->Top.Y == Y && e->NextInLML;
	}
	//------------------------------------------------------------------------------
	//------------------------------------------------------------------------------
	//------------------------------------------------------------------------------

	void ClipperEx::SwapPositionsInSEL(TEdge* Edge1, TEdge* Edge2)
	{
		if (!(Edge1->NextInSEL) && !(Edge1->PrevInSEL)) return;
		if (!(Edge2->NextInSEL) && !(Edge2->PrevInSEL)) return;

		if (Edge1->NextInSEL == Edge2)
		{
			TEdge* Next = Edge2->NextInSEL;
			if (Next) Next->PrevInSEL = Edge1;
			TEdge* Prev = Edge1->PrevInSEL;
			if (Prev) Prev->NextInSEL = Edge2;
			Edge2->PrevInSEL = Prev;
			Edge2->NextInSEL = Edge1;
			Edge1->PrevInSEL = Edge2;
			Edge1->NextInSEL = Next;
		}
		else if (Edge2->NextInSEL == Edge1)
		{
			TEdge* Next = Edge1->NextInSEL;
			if (Next) Next->PrevInSEL = Edge2;
			TEdge* Prev = Edge2->PrevInSEL;
			if (Prev) Prev->NextInSEL = Edge1;
			Edge1->PrevInSEL = Prev;
			Edge1->NextInSEL = Edge2;
			Edge2->PrevInSEL = Edge1;
			Edge2->NextInSEL = Next;
		}
		else
		{
			TEdge* Next = Edge1->NextInSEL;
			TEdge* Prev = Edge1->PrevInSEL;
			Edge1->NextInSEL = Edge2->NextInSEL;
			if (Edge1->NextInSEL) Edge1->NextInSEL->PrevInSEL = Edge1;
			Edge1->PrevInSEL = Edge2->PrevInSEL;
			if (Edge1->PrevInSEL) Edge1->PrevInSEL->NextInSEL = Edge1;
			Edge2->NextInSEL = Next;
			if (Edge2->NextInSEL) Edge2->NextInSEL->PrevInSEL = Edge2;
			Edge2->PrevInSEL = Prev;
			if (Edge2->PrevInSEL) Edge2->PrevInSEL->NextInSEL = Edge2;
		}

		if (!Edge1->PrevInSEL) m_SortedEdges = Edge1;
		else if (!Edge2->PrevInSEL) m_SortedEdges = Edge2;
	}
	//------------------------------------------------------------------------------
	//------------------------------------------------------------------------

	/*******************************************************************************
	* Notes: Horizontal edges (HEs) at scanline intersections (ie at the Top or    *
	* Bottom of a scanbeam) are processed as if layered. The order in which HEs    *
	* are processed doesn't matter. HEs intersect with other HE Bot.Xs only [#]    *
	* (or they could intersect with Top.Xs only, ie EITHER Bot.Xs OR Top.Xs),      *
	* and with other non-horizontal edges [*]. Once these intersections are        *
	* processed, intermediate HEs then 'promote' the Edge above (NextInLML) into   *
	* the AEL. These 'promoted' edges may in turn intersect [%] with other HEs.    *
	*******************************************************************************/

	void ClipperEx::ProcessHorizontal(TEdge* horzEdge)
	{
		Direction dir;
		cInt horzLeft, horzRight;
		bool IsOpen = (horzEdge->WindDelta == 0);

		GetHorzDirection(*horzEdge, dir, horzLeft, horzRight);

		TEdge* eLastHorz = horzEdge, * eMaxPair = 0;
		while (eLastHorz->NextInLML && IsHorizontal(*eLastHorz->NextInLML))
			eLastHorz = eLastHorz->NextInLML;
		if (!eLastHorz->NextInLML)
			eMaxPair = GetMaximaPair(eLastHorz);

		MaximaList::const_iterator maxIt;
		MaximaList::const_reverse_iterator maxRit;
		if (m_Maxima.size() > 0)
		{
			//get the first maxima in range (X) ...
			if (dir == dLeftToRight)
			{
				maxIt = m_Maxima.begin();
				while (maxIt != m_Maxima.end() && *maxIt <= horzEdge->Bot.X) maxIt++;
				if (maxIt != m_Maxima.end() && *maxIt >= eLastHorz->Top.X)
					maxIt = m_Maxima.end();
			}
			else
			{
				maxRit = m_Maxima.rbegin();
				while (maxRit != m_Maxima.rend() && *maxRit > horzEdge->Bot.X) maxRit++;
				if (maxRit != m_Maxima.rend() && *maxRit <= eLastHorz->Top.X)
					maxRit = m_Maxima.rend();
			}
		}

		OutPt* op1 = 0;

		for (;;) //loop through consec. horizontal edges
		{

			bool IsLastHorz = (horzEdge == eLastHorz);
			TEdge* e = GetNextInAEL(horzEdge, dir);
			while (e)
			{

				//this code block inserts extra coords into horizontal edges (in output
				//polygons) whereever maxima touch these horizontal edges. This helps
				//'simplifying' polygons (ie if the Simplify property is set).
				if (m_Maxima.size() > 0)
				{
					if (dir == dLeftToRight)
					{
						while (maxIt != m_Maxima.end() && *maxIt < e->Curr.X)
						{
							if (horzEdge->OutIdx >= 0 && !IsOpen)
								AddOutPt(horzEdge, IntPoint(*maxIt, horzEdge->Bot.Y));
							maxIt++;
						}
					}
					else
					{
						while (maxRit != m_Maxima.rend() && *maxRit > e->Curr.X)
						{
							if (horzEdge->OutIdx >= 0 && !IsOpen)
								AddOutPt(horzEdge, IntPoint(*maxRit, horzEdge->Bot.Y));
							maxRit++;
						}
					}
				};

				if ((dir == dLeftToRight && e->Curr.X > horzRight) ||
					(dir == dRightToLeft && e->Curr.X < horzLeft)) break;

				//Also break if we've got to the end of an intermediate horizontal edge ...
				//nb: Smaller Dx's are to the right of larger Dx's ABOVE the horizontal.
				if (e->Curr.X == horzEdge->Top.X && horzEdge->NextInLML &&
					e->Dx < horzEdge->NextInLML->Dx) break;

				if (horzEdge->OutIdx >= 0 && !IsOpen)  //note: may be done multiple times
				{
#ifdef use_xyz
					if (dir == dLeftToRight) SetZ(e->Curr, *horzEdge, *e);
					else SetZ(e->Curr, *e, *horzEdge);
#endif      
					op1 = AddOutPt(horzEdge, e->Curr);
					TEdge* eNextHorz = m_SortedEdges;
					while (eNextHorz)
					{
						if (eNextHorz->OutIdx >= 0 &&
							HorzSegmentsOverlap(horzEdge->Bot.X,
								horzEdge->Top.X, eNextHorz->Bot.X, eNextHorz->Top.X))
						{
							OutPt* op2 = GetLastOutPt(eNextHorz);
							AddJoin(op2, op1, eNextHorz->Top);
						}
						eNextHorz = eNextHorz->NextInSEL;
					}
					AddGhostJoin(op1, horzEdge->Bot);
				}

				//OK, so far we're still in range of the horizontal Edge  but make sure
				//we're at the last of consec. horizontals when matching with eMaxPair
				if (e == eMaxPair && IsLastHorz)
				{
					if (horzEdge->OutIdx >= 0)
						AddLocalMaxPoly(horzEdge, eMaxPair, horzEdge->Top);
					DeleteFromAEL(horzEdge);
					DeleteFromAEL(eMaxPair);
					return;
				}

				if (dir == dLeftToRight)
				{
					IntPoint Pt = IntPoint(e->Curr.X, horzEdge->Curr.Y);
					IntersectEdges(horzEdge, e, Pt);
				}
				else
				{
					IntPoint Pt = IntPoint(e->Curr.X, horzEdge->Curr.Y);
					IntersectEdges(e, horzEdge, Pt);
				}
				TEdge* eNext = GetNextInAEL(e, dir);
				SwapPositionsInAEL(horzEdge, e);
				e = eNext;
			} //end while(e)

			//Break out of loop if HorzEdge.NextInLML is not also horizontal ...
			if (!horzEdge->NextInLML || !IsHorizontal(*horzEdge->NextInLML)) break;

			UpdateEdgeIntoAEL(horzEdge);
			if (horzEdge->OutIdx >= 0) AddOutPt(horzEdge, horzEdge->Bot);
			GetHorzDirection(*horzEdge, dir, horzLeft, horzRight);

		} //end for (;;)

		if (horzEdge->OutIdx >= 0 && !op1)
		{
			op1 = GetLastOutPt(horzEdge);
			TEdge* eNextHorz = m_SortedEdges;
			while (eNextHorz)
			{
				if (eNextHorz->OutIdx >= 0 &&
					HorzSegmentsOverlap(horzEdge->Bot.X,
						horzEdge->Top.X, eNextHorz->Bot.X, eNextHorz->Top.X))
				{
					OutPt* op2 = GetLastOutPt(eNextHorz);
					AddJoin(op2, op1, eNextHorz->Top);
				}
				eNextHorz = eNextHorz->NextInSEL;
			}
			AddGhostJoin(op1, horzEdge->Top);
		}

		if (horzEdge->NextInLML)
		{
			if (horzEdge->OutIdx >= 0)
			{
				op1 = AddOutPt(horzEdge, horzEdge->Top);
				UpdateEdgeIntoAEL(horzEdge);
				if (horzEdge->WindDelta == 0) return;
				//nb: HorzEdge is no longer horizontal here
				TEdge* ePrev = horzEdge->PrevInAEL;
				TEdge* eNext = horzEdge->NextInAEL;
				if (ePrev && ePrev->Curr.X == horzEdge->Bot.X &&
					ePrev->Curr.Y == horzEdge->Bot.Y && ePrev->WindDelta != 0 &&
					(ePrev->OutIdx >= 0 && ePrev->Curr.Y > ePrev->Top.Y &&
						SlopesEqual(*horzEdge, *ePrev, m_UseFullRange)))
				{
					OutPt* op2 = AddOutPt(ePrev, horzEdge->Bot);
					AddJoin(op1, op2, horzEdge->Top);
				}
				else if (eNext && eNext->Curr.X == horzEdge->Bot.X &&
					eNext->Curr.Y == horzEdge->Bot.Y && eNext->WindDelta != 0 &&
					eNext->OutIdx >= 0 && eNext->Curr.Y > eNext->Top.Y &&
					SlopesEqual(*horzEdge, *eNext, m_UseFullRange))
				{
					OutPt* op2 = AddOutPt(eNext, horzEdge->Bot);
					AddJoin(op1, op2, horzEdge->Top);
				}
			}
			else
				UpdateEdgeIntoAEL(horzEdge);
		}
		else
		{
			if (horzEdge->OutIdx >= 0) AddOutPt(horzEdge, horzEdge->Top);
			DeleteFromAEL(horzEdge);
		}
	}
	//------------------------------------------------------------------------------

	bool ClipperEx::ProcessIntersections(const cInt topY)
	{
		if (!m_ActiveEdges) return true;
		try {
			BuildIntersectList(topY);
			size_t IlSize = m_IntersectList.size();
			if (IlSize == 0) return true;
			if (IlSize == 1 || FixupIntersectionOrder()) ProcessIntersectList();
			else return false;
		}
		catch (...)
		{
			m_SortedEdges = 0;
			DisposeIntersectNodes();
			throw clipperException("ProcessIntersections error");
		}
		m_SortedEdges = 0;
		return true;
	}
	//------------------------------------------------------------------------------

	void ClipperEx::DisposeIntersectNodes()
	{
		for (size_t i = 0; i < m_IntersectList.size(); ++i)
			delete m_IntersectList[i];
		m_IntersectList.clear();
	}
	//------------------------------------------------------------------------------

	void ClipperEx::BuildIntersectList(const cInt topY)
	{
		if (!m_ActiveEdges) return;

		//prepare for sorting ...
		TEdge* e = m_ActiveEdges;
		m_SortedEdges = e;
		while (e)
		{
			e->PrevInSEL = e->PrevInAEL;
			e->NextInSEL = e->NextInAEL;
			e->Curr.X = TopX(*e, topY);
			e = e->NextInAEL;
		}

		//bubblesort ...
		bool isModified;
		do
		{
			isModified = false;
			e = m_SortedEdges;
			while (e->NextInSEL)
			{
				TEdge* eNext = e->NextInSEL;
				IntPoint Pt;
				if (e->Curr.X > eNext->Curr.X)
				{
					IntersectPoint(*e, *eNext, Pt);
					if (Pt.Y < topY) Pt = IntPoint(TopX(*e, topY), topY);
					IntersectNode* newNode = new IntersectNode;
					newNode->Edge1 = e;
					newNode->Edge2 = eNext;
					newNode->Pt = Pt;
					m_IntersectList.push_back(newNode);

					SwapPositionsInSEL(e, eNext);
					isModified = true;
				}
				else
					e = eNext;
			}
			if (e->PrevInSEL) e->PrevInSEL->NextInSEL = 0;
			else break;
		} while (isModified);
		m_SortedEdges = 0; //important
	}
	//------------------------------------------------------------------------------


	void ClipperEx::ProcessIntersectList()
	{
		for (size_t i = 0; i < m_IntersectList.size(); ++i)
		{
			IntersectNode* iNode = m_IntersectList[i];
			{
				IntersectEdges(iNode->Edge1, iNode->Edge2, iNode->Pt);
				SwapPositionsInAEL(iNode->Edge1, iNode->Edge2);
			}
			delete iNode;
		}
		m_IntersectList.clear();
	}
	//------------------------------------------------------------------------------
	//------------------------------------------------------------------------------

	inline bool EdgesAdjacent(const IntersectNode& inode)
	{
		return (inode.Edge1->NextInSEL == inode.Edge2) ||
			(inode.Edge1->PrevInSEL == inode.Edge2);
	}
	//------------------------------------------------------------------------------

	bool ClipperEx::FixupIntersectionOrder()
	{
		//pre-condition: intersections are sorted Bottom-most first.
		//Now it's crucial that intersections are made only between adjacent edges,
		//so to ensure this the order of intersections may need adjusting ...
		CopyAELToSEL();
		std::sort(m_IntersectList.begin(), m_IntersectList.end(), IntersectListSort);
		size_t cnt = m_IntersectList.size();
		for (size_t i = 0; i < cnt; ++i)
		{
			if (!EdgesAdjacent(*m_IntersectList[i]))
			{
				size_t j = i + 1;
				while (j < cnt && !EdgesAdjacent(*m_IntersectList[j])) j++;
				if (j == cnt)  return false;
				std::swap(m_IntersectList[i], m_IntersectList[j]);
			}
			SwapPositionsInSEL(m_IntersectList[i]->Edge1, m_IntersectList[i]->Edge2);
		}
		return true;
	}
	//------------------------------------------------------------------------------

	void ClipperEx::DoMaxima(TEdge* e)
	{
		TEdge* eMaxPair = GetMaximaPairEx(e);
		if (!eMaxPair)
		{
			if (e->OutIdx >= 0)
				AddOutPt(e, e->Top);
			DeleteFromAEL(e);
			return;
		}

		TEdge* eNext = e->NextInAEL;
		while (eNext && eNext != eMaxPair)
		{
			IntersectEdges(e, eNext, e->Top);
			SwapPositionsInAEL(e, eNext);
			eNext = e->NextInAEL;
		}

		if (e->OutIdx == Unassigned && eMaxPair->OutIdx == Unassigned)
		{
			DeleteFromAEL(e);
			DeleteFromAEL(eMaxPair);
		}
		else if (e->OutIdx >= 0 && eMaxPair->OutIdx >= 0)
		{
			if (e->OutIdx >= 0) AddLocalMaxPoly(e, eMaxPair, e->Top);
			DeleteFromAEL(e);
			DeleteFromAEL(eMaxPair);
		}
#ifdef use_lines
		else if (e->WindDelta == 0)
		{
			if (e->OutIdx >= 0)
			{
				AddOutPt(e, e->Top);
				e->OutIdx = Unassigned;
			}
			DeleteFromAEL(e);

			if (eMaxPair->OutIdx >= 0)
			{
				AddOutPt(eMaxPair, e->Top);
				eMaxPair->OutIdx = Unassigned;
			}
			DeleteFromAEL(eMaxPair);
		}
#endif
		else throw clipperException("DoMaxima error");
	}
	//------------------------------------------------------------------------------

	void ClipperEx::ProcessEdgesAtTopOfScanbeam(const cInt topY)
	{
		TEdge* e = m_ActiveEdges;
		while (e)
		{
			//1. process maxima, treating them as if they're 'bent' horizontal edges,
			//   but exclude maxima with horizontal edges. nb: e can't be a horizontal.
			bool IsMaximaEdge = IsMaxima(e, topY);

			if (IsMaximaEdge)
			{
				TEdge* eMaxPair = GetMaximaPairEx(e);
				IsMaximaEdge = (!eMaxPair || !IsHorizontal(*eMaxPair));
			}

			if (IsMaximaEdge)
			{
				if (m_StrictSimple) m_Maxima.push_back(e->Top.X);
				TEdge* ePrev = e->PrevInAEL;
				DoMaxima(e);
				if (!ePrev) e = m_ActiveEdges;
				else e = ePrev->NextInAEL;
			}
			else
			{
				//2. promote horizontal edges, otherwise update Curr.X and Curr.Y ...
				if (IsIntermediate(e, topY) && IsHorizontal(*e->NextInLML))
				{
					UpdateEdgeIntoAEL(e);
					if (e->OutIdx >= 0)
						AddOutPt(e, e->Bot);
					AddEdgeToSEL(e);
				}
				else
				{
					e->Curr.X = TopX(*e, topY);
					e->Curr.Y = topY;
#ifdef use_xyz
					e->Curr.Z = topY == e->Top.Y ? e->Top.Z : (topY == e->Bot.Y ? e->Bot.Z : 0);
#endif
				}

				//When StrictlySimple and 'e' is being touched by another edge, then
				//make sure both edges have a vertex here ...
				if (m_StrictSimple)
				{
					TEdge* ePrev = e->PrevInAEL;
					if ((e->OutIdx >= 0) && (e->WindDelta != 0) && ePrev && (ePrev->OutIdx >= 0) &&
						(ePrev->Curr.X == e->Curr.X) && (ePrev->WindDelta != 0))
					{
						IntPoint pt = e->Curr;
#ifdef use_xyz
						SetZ(pt, *ePrev, *e);
#endif
						OutPt* op = AddOutPt(ePrev, pt);
						OutPt* op2 = AddOutPt(e, pt);
						AddJoin(op, op2, pt); //StrictlySimple (type-3) join
					}
				}

				e = e->NextInAEL;
			}
		}

		//3. Process horizontals at the Top of the scanbeam ...
		m_Maxima.sort();
		ProcessHorizontals();
		m_Maxima.clear();

		//4. Promote intermediate vertices ...
		e = m_ActiveEdges;
		while (e)
		{
			if (IsIntermediate(e, topY))
			{
				OutPt* op = 0;
				if (e->OutIdx >= 0)
					op = AddOutPt(e, e->Top);
				UpdateEdgeIntoAEL(e);

				//if output polygons share an edge, they'll need joining later ...
				TEdge* ePrev = e->PrevInAEL;
				TEdge* eNext = e->NextInAEL;
				if (ePrev && ePrev->Curr.X == e->Bot.X &&
					ePrev->Curr.Y == e->Bot.Y && op &&
					ePrev->OutIdx >= 0 && ePrev->Curr.Y > ePrev->Top.Y &&
					SlopesEqual(e->Curr, e->Top, ePrev->Curr, ePrev->Top, m_UseFullRange) &&
					(e->WindDelta != 0) && (ePrev->WindDelta != 0))
				{
					OutPt* op2 = AddOutPt(ePrev, e->Bot);
					AddJoin(op, op2, e->Top);
				}
				else if (eNext && eNext->Curr.X == e->Bot.X &&
					eNext->Curr.Y == e->Bot.Y && op &&
					eNext->OutIdx >= 0 && eNext->Curr.Y > eNext->Top.Y &&
					SlopesEqual(e->Curr, e->Top, eNext->Curr, eNext->Top, m_UseFullRange) &&
					(e->WindDelta != 0) && (eNext->WindDelta != 0))
				{
					OutPt* op2 = AddOutPt(eNext, e->Bot);
					AddJoin(op, op2, e->Top);
				}
			}
			e = e->NextInAEL;
		}
	}
	//------------------------------------------------------------------------------

	void ClipperEx::FixupOutPolyline(OutRec& outrec)
	{
		OutPt* pp = outrec.Pts;
		OutPt* lastPP = pp->Prev;
		while (pp != lastPP)
		{
			pp = pp->Next;
			if (pp->Pt == pp->Prev->Pt)
			{
				if (pp == lastPP) lastPP = pp->Prev;
				OutPt* tmpPP = pp->Prev;
				tmpPP->Next = pp->Next;
				pp->Next->Prev = tmpPP;
				delete pp;
				pp = tmpPP;
			}
		}

		if (pp == pp->Prev)
		{
			DisposeOutPts(pp);
			outrec.Pts = 0;
			return;
		}
	}
	//------------------------------------------------------------------------------

	void ClipperEx::FixupOutPolygon(OutRec& outrec)
	{
		//FixupOutPolygon() - removes duplicate points and simplifies consecutive
		//parallel edges by removing the middle vertex.
		OutPt* lastOK = 0;
		outrec.BottomPt = 0;
		OutPt* pp = outrec.Pts;
		bool preserveCol = m_PreserveCollinear || m_StrictSimple;

		for (;;)
		{
			if (pp->Prev == pp || pp->Prev == pp->Next)
			{
				DisposeOutPts(pp);
				outrec.Pts = 0;
				return;
			}

			//test for duplicate points and collinear edges ...
			if ((pp->Pt == pp->Next->Pt) || (pp->Pt == pp->Prev->Pt) ||
				(SlopesEqual(pp->Prev->Pt, pp->Pt, pp->Next->Pt, m_UseFullRange) &&
				(!preserveCol || !Pt2IsBetweenPt1AndPt3(pp->Prev->Pt, pp->Pt, pp->Next->Pt))))
			{
				lastOK = 0;
				OutPt* tmp = pp;
				pp->Prev->Next = pp->Next;
				pp->Next->Prev = pp->Prev;
				pp = pp->Prev;
				delete tmp;
			}
			else if (pp == lastOK) break;
			else
			{
				if (!lastOK) lastOK = pp;
				pp = pp->Next;
			}
		}
		outrec.Pts = pp;
	}
	//------------------------------------------------------------------------------
	//------------------------------------------------------------------------------

	void ClipperEx::BuildResult(Paths& polys)
	{
		polys.reserve(m_PolyOuts.size());
		for (PolyOutList::size_type i = 0; i < m_PolyOuts.size(); ++i)
		{
			if (!m_PolyOuts[i]->Pts) continue;
			Path pg;
			OutPt* p = m_PolyOuts[i]->Pts->Prev;
			int cnt = PointCount(p);
			if (cnt < 2) continue;
			pg.reserve(cnt);
			for (int i = 0; i < cnt; ++i)
			{
				pg.push_back(p->Pt);
				p = p->Prev;
			}
			polys.push_back(pg);
		}
	}
	//------------------------------------------------------------------------------

	void ClipperEx::BuildResult2(PolyTree& polytree)
	{
		polytree.Clear();
		polytree.AllNodes.reserve(m_PolyOuts.size());
		//add each output polygon/contour to polytree ...
		for (PolyOutList::size_type i = 0; i < m_PolyOuts.size(); i++)
		{
			OutRec* outRec = m_PolyOuts[i];
			int cnt = PointCount(outRec->Pts);
			if ((outRec->IsOpen && cnt < 2) || (!outRec->IsOpen && cnt < 3)) continue;
			FixHoleLinkage(*outRec);
			PolyNode* pn = new PolyNode();
			//nb: polytree takes ownership of all the PolyNodes
			polytree.AllNodes.push_back(pn);
			outRec->PolyNd = pn;
			pn->Parent = 0;
			pn->Index = 0;
			pn->Contour.reserve(cnt);
			OutPt* op = outRec->Pts->Prev;
			for (int j = 0; j < cnt; j++)
			{
				pn->Contour.push_back(op->Pt);
				op = op->Prev;
			}
		}

		//fixup PolyNode links etc ...
		polytree.Childs.reserve(m_PolyOuts.size());
		for (PolyOutList::size_type i = 0; i < m_PolyOuts.size(); i++)
		{
			OutRec* outRec = m_PolyOuts[i];
			if (!outRec->PolyNd) continue;
			if (outRec->IsOpen)
			{
				outRec->PolyNd->m_IsOpen = true;
				polytree.AddChild(*outRec->PolyNd);
			}
			else if (outRec->FirstLeft && outRec->FirstLeft->PolyNd)
				outRec->FirstLeft->PolyNd->AddChild(*outRec->PolyNd);
			else
				polytree.AddChild(*outRec->PolyNd);
		}
	}
	//------------------------------------------------------------------------------
	//------------------------------------------------------------------------------

	inline bool E2InsertsBeforeE1(TEdge& e1, TEdge& e2)
	{
		if (e2.Curr.X == e1.Curr.X)
		{
			if (e2.Top.Y > e1.Top.Y)
				return e2.Top.X < TopX(e1, e2.Top.Y);
			else return e1.Top.X > TopX(e2, e1.Top.Y);
		}
		else return e2.Curr.X < e1.Curr.X;
	}
	//------------------------------------------------------------------------------

	//------------------------------------------------------------------------------

	inline void UpdateOutPtIdxs(OutRec& outrec)
	{
		OutPt* op = outrec.Pts;
		do
		{
			op->Idx = outrec.Idx;
			op = op->Prev;
		} while (op != outrec.Pts);
	}
	//------------------------------------------------------------------------------

	void ClipperEx::InsertEdgeIntoAEL(TEdge* edge, TEdge* startEdge)
	{
		if (!m_ActiveEdges)
		{
			edge->PrevInAEL = 0;
			edge->NextInAEL = 0;
			m_ActiveEdges = edge;
		}
		else if (!startEdge && E2InsertsBeforeE1(*m_ActiveEdges, *edge))
		{
			edge->PrevInAEL = 0;
			edge->NextInAEL = m_ActiveEdges;
			m_ActiveEdges->PrevInAEL = edge;
			m_ActiveEdges = edge;
		}
		else
		{
			if (!startEdge) startEdge = m_ActiveEdges;
			while (startEdge->NextInAEL &&
				!E2InsertsBeforeE1(*startEdge->NextInAEL, *edge))
				startEdge = startEdge->NextInAEL;
			edge->NextInAEL = startEdge->NextInAEL;
			if (startEdge->NextInAEL) startEdge->NextInAEL->PrevInAEL = edge;
			edge->PrevInAEL = startEdge;
			startEdge->NextInAEL = edge;
		}
	}
	//----------------------------------------------------------------------
	void ClipperEx::DoSimplePolygons()
	{
		PolyOutList::size_type i = 0;
		while (i < m_PolyOuts.size())
		{
			OutRec* outrec = m_PolyOuts[i++];
			OutPt* op = outrec->Pts;
			if (!op || outrec->IsOpen) continue;
			do //for each Pt in Polygon until duplicate found do ...
			{
				OutPt* op2 = op->Next;
				while (op2 != outrec->Pts)
				{
					if ((op->Pt == op2->Pt) && op2->Next != op && op2->Prev != op)
					{
						//split the polygon into two ...
						OutPt* op3 = op->Prev;
						OutPt* op4 = op2->Prev;
						op->Prev = op4;
						op4->Next = op;
						op2->Prev = op3;
						op3->Next = op2;

						outrec->Pts = op;
						OutRec* outrec2 = CreateOutRec();
						outrec2->Pts = op2;
						UpdateOutPtIdxs(*outrec2);
						if (Poly2ContainsPoly1(outrec2->Pts, outrec->Pts))
						{
							//OutRec2 is contained by OutRec1 ...
							outrec2->IsHole = !outrec->IsHole;
							outrec2->FirstLeft = outrec;
							if (m_UsingPolyTree) FixupFirstLefts2(outrec2, outrec);
						}
						else
							if (Poly2ContainsPoly1(outrec->Pts, outrec2->Pts))
							{
								//OutRec1 is contained by OutRec2 ...
								outrec2->IsHole = outrec->IsHole;
								outrec->IsHole = !outrec2->IsHole;
								outrec2->FirstLeft = outrec->FirstLeft;
								outrec->FirstLeft = outrec2;
								if (m_UsingPolyTree) FixupFirstLefts2(outrec, outrec2);
							}
							else
							{
								//the 2 polygons are separate ...
								outrec2->IsHole = outrec->IsHole;
								outrec2->FirstLeft = outrec->FirstLeft;
								if (m_UsingPolyTree) FixupFirstLefts1(outrec, outrec2);
							}
						op2 = op; //ie get ready for the Next iteration
					}
					op2 = op2->Next;
				}
				op = op->Next;
			} while (op != outrec->Pts);
		}
	}
	//------------------------------------------------------------------------------

	bool ClipperEx::JoinPoints(Join* j, OutRec* outRec1, OutRec* outRec2)
	{
		OutPt* op1 = j->OutPt1, * op1b;
		OutPt* op2 = j->OutPt2, * op2b;

		//There are 3 kinds of joins for output polygons ...
		//1. Horizontal joins where Join.OutPt1 & Join.OutPt2 are vertices anywhere
		//along (horizontal) collinear edges (& Join.OffPt is on the same horizontal).
		//2. Non-horizontal joins where Join.OutPt1 & Join.OutPt2 are at the same
		//location at the Bottom of the overlapping segment (& Join.OffPt is above).
		//3. StrictSimple joins where edges touch but are not collinear and where
		//Join.OutPt1, Join.OutPt2 & Join.OffPt all share the same point.
		bool isHorizontal = (j->OutPt1->Pt.Y == j->OffPt.Y);

		if (isHorizontal && (j->OffPt == j->OutPt1->Pt) &&
			(j->OffPt == j->OutPt2->Pt))
		{
			//Strictly Simple join ...
			if (outRec1 != outRec2) return false;
			op1b = j->OutPt1->Next;
			while (op1b != op1 && (op1b->Pt == j->OffPt))
				op1b = op1b->Next;
			bool reverse1 = (op1b->Pt.Y > j->OffPt.Y);
			op2b = j->OutPt2->Next;
			while (op2b != op2 && (op2b->Pt == j->OffPt))
				op2b = op2b->Next;
			bool reverse2 = (op2b->Pt.Y > j->OffPt.Y);
			if (reverse1 == reverse2) return false;
			if (reverse1)
			{
				op1b = DupOutPt(op1, false);
				op2b = DupOutPt(op2, true);
				op1->Prev = op2;
				op2->Next = op1;
				op1b->Next = op2b;
				op2b->Prev = op1b;
				j->OutPt1 = op1;
				j->OutPt2 = op1b;
				return true;
			}
			else
			{
				op1b = DupOutPt(op1, true);
				op2b = DupOutPt(op2, false);
				op1->Next = op2;
				op2->Prev = op1;
				op1b->Prev = op2b;
				op2b->Next = op1b;
				j->OutPt1 = op1;
				j->OutPt2 = op1b;
				return true;
			}
		}
		else if (isHorizontal)
		{
			//treat horizontal joins differently to non-horizontal joins since with
			//them we're not yet sure where the overlapping is. OutPt1.Pt & OutPt2.Pt
			//may be anywhere along the horizontal edge.
			op1b = op1;
			while (op1->Prev->Pt.Y == op1->Pt.Y && op1->Prev != op1b && op1->Prev != op2)
				op1 = op1->Prev;
			while (op1b->Next->Pt.Y == op1b->Pt.Y && op1b->Next != op1 && op1b->Next != op2)
				op1b = op1b->Next;
			if (op1b->Next == op1 || op1b->Next == op2) return false; //a flat 'polygon'

			op2b = op2;
			while (op2->Prev->Pt.Y == op2->Pt.Y && op2->Prev != op2b && op2->Prev != op1b)
				op2 = op2->Prev;
			while (op2b->Next->Pt.Y == op2b->Pt.Y && op2b->Next != op2 && op2b->Next != op1)
				op2b = op2b->Next;
			if (op2b->Next == op2 || op2b->Next == op1) return false; //a flat 'polygon'

			cInt Left, Right;
			//Op1 --> Op1b & Op2 --> Op2b are the extremites of the horizontal edges
			if (!GetOverlap(op1->Pt.X, op1b->Pt.X, op2->Pt.X, op2b->Pt.X, Left, Right))
				return false;

			//DiscardLeftSide: when overlapping edges are joined, a spike will created
			//which needs to be cleaned up. However, we don't want Op1 or Op2 caught up
			//on the discard Side as either may still be needed for other joins ...
			IntPoint Pt;
			bool DiscardLeftSide;
			if (op1->Pt.X >= Left && op1->Pt.X <= Right)
			{
				Pt = op1->Pt; DiscardLeftSide = (op1->Pt.X > op1b->Pt.X);
			}
			else if (op2->Pt.X >= Left && op2->Pt.X <= Right)
			{
				Pt = op2->Pt; DiscardLeftSide = (op2->Pt.X > op2b->Pt.X);
			}
			else if (op1b->Pt.X >= Left && op1b->Pt.X <= Right)
			{
				Pt = op1b->Pt; DiscardLeftSide = op1b->Pt.X > op1->Pt.X;
			}
			else
			{
				Pt = op2b->Pt; DiscardLeftSide = (op2b->Pt.X > op2->Pt.X);
			}
			j->OutPt1 = op1; j->OutPt2 = op2;
			return JoinHorz(op1, op1b, op2, op2b, Pt, DiscardLeftSide);
		}
		else
		{
			//nb: For non-horizontal joins ...
			//    1. Jr.OutPt1.Pt.Y == Jr.OutPt2.Pt.Y
			//    2. Jr.OutPt1.Pt > Jr.OffPt.Y

			//make sure the polygons are correctly oriented ...
			op1b = op1->Next;
			while ((op1b->Pt == op1->Pt) && (op1b != op1)) op1b = op1b->Next;
			bool Reverse1 = ((op1b->Pt.Y > op1->Pt.Y) ||
				!SlopesEqual(op1->Pt, op1b->Pt, j->OffPt, m_UseFullRange));
			if (Reverse1)
			{
				op1b = op1->Prev;
				while ((op1b->Pt == op1->Pt) && (op1b != op1)) op1b = op1b->Prev;
				if ((op1b->Pt.Y > op1->Pt.Y) ||
					!SlopesEqual(op1->Pt, op1b->Pt, j->OffPt, m_UseFullRange)) return false;
			};
			op2b = op2->Next;
			while ((op2b->Pt == op2->Pt) && (op2b != op2))op2b = op2b->Next;
			bool Reverse2 = ((op2b->Pt.Y > op2->Pt.Y) ||
				!SlopesEqual(op2->Pt, op2b->Pt, j->OffPt, m_UseFullRange));
			if (Reverse2)
			{
				op2b = op2->Prev;
				while ((op2b->Pt == op2->Pt) && (op2b != op2)) op2b = op2b->Prev;
				if ((op2b->Pt.Y > op2->Pt.Y) ||
					!SlopesEqual(op2->Pt, op2b->Pt, j->OffPt, m_UseFullRange)) return false;
			}

			if ((op1b == op1) || (op2b == op2) || (op1b == op2b) ||
				((outRec1 == outRec2) && (Reverse1 == Reverse2))) return false;

			if (Reverse1)
			{
				op1b = DupOutPt(op1, false);
				op2b = DupOutPt(op2, true);
				op1->Prev = op2;
				op2->Next = op1;
				op1b->Next = op2b;
				op2b->Prev = op1b;
				j->OutPt1 = op1;
				j->OutPt2 = op1b;
				return true;
			}
			else
			{
				op1b = DupOutPt(op1, true);
				op2b = DupOutPt(op2, false);
				op1->Next = op2;
				op2->Prev = op1;
				op1b->Prev = op2b;
				op2b->Next = op1b;
				j->OutPt1 = op1;
				j->OutPt2 = op1b;
				return true;
			}
		}
	}
	//----------------------------------------------------------------------

	static OutRec* ParseFirstLeft(OutRec* FirstLeft)
	{
		while (FirstLeft && !FirstLeft->Pts)
			FirstLeft = FirstLeft->FirstLeft;
		return FirstLeft;
	}
	//------------------------------------------------------------------------------

	void ClipperEx::FixupFirstLefts1(OutRec* OldOutRec, OutRec* NewOutRec)
	{
		//tests if NewOutRec contains the polygon before reassigning FirstLeft
		for (PolyOutList::size_type i = 0; i < m_PolyOuts.size(); ++i)
		{
			OutRec* outRec = m_PolyOuts[i];
			OutRec* firstLeft = ParseFirstLeft(outRec->FirstLeft);
			if (outRec->Pts && firstLeft == OldOutRec)
			{
				if (Poly2ContainsPoly1(outRec->Pts, NewOutRec->Pts))
					outRec->FirstLeft = NewOutRec;
			}
		}
	}
	//----------------------------------------------------------------------

	void ClipperEx::FixupFirstLefts2(OutRec* InnerOutRec, OutRec* OuterOutRec)
	{
		//A polygon has split into two such that one is now the inner of the other.
		//It's possible that these polygons now wrap around other polygons, so check
		//every polygon that's also contained by OuterOutRec's FirstLeft container
		//(including 0) to see if they've become inner to the new inner polygon ...
		OutRec* orfl = OuterOutRec->FirstLeft;
		for (PolyOutList::size_type i = 0; i < m_PolyOuts.size(); ++i)
		{
			OutRec* outRec = m_PolyOuts[i];

			if (!outRec->Pts || outRec == OuterOutRec || outRec == InnerOutRec)
				continue;
			OutRec* firstLeft = ParseFirstLeft(outRec->FirstLeft);
			if (firstLeft != orfl && firstLeft != InnerOutRec && firstLeft != OuterOutRec)
				continue;
			if (Poly2ContainsPoly1(outRec->Pts, InnerOutRec->Pts))
				outRec->FirstLeft = InnerOutRec;
			else if (Poly2ContainsPoly1(outRec->Pts, OuterOutRec->Pts))
				outRec->FirstLeft = OuterOutRec;
			else if (outRec->FirstLeft == InnerOutRec || outRec->FirstLeft == OuterOutRec)
				outRec->FirstLeft = orfl;
		}
	}
	//----------------------------------------------------------------------
	void ClipperEx::FixupFirstLefts3(OutRec* OldOutRec, OutRec* NewOutRec)
	{
		//reassigns FirstLeft WITHOUT testing if NewOutRec contains the polygon
		for (PolyOutList::size_type i = 0; i < m_PolyOuts.size(); ++i)
		{
			OutRec* outRec = m_PolyOuts[i];
			OutRec* firstLeft = ParseFirstLeft(outRec->FirstLeft);
			if (outRec->Pts && firstLeft == OldOutRec)
				outRec->FirstLeft = NewOutRec;
		}
	}
	//----------------------------------------------------------------------

	void ClipperEx::JoinCommonEdges()
	{
		for (JoinList::size_type i = 0; i < m_Joins.size(); i++)
		{
			Join* join = m_Joins[i];

			OutRec* outRec1 = GetOutRec(join->OutPt1->Idx);
			OutRec* outRec2 = GetOutRec(join->OutPt2->Idx);

			if (!outRec1->Pts || !outRec2->Pts) continue;
			if (outRec1->IsOpen || outRec2->IsOpen) continue;

			//get the polygon fragment with the correct hole state (FirstLeft)
			//before calling JoinPoints() ...
			OutRec* holeStateRec;
			if (outRec1 == outRec2) holeStateRec = outRec1;
			else if (OutRec1RightOfOutRec2(outRec1, outRec2)) holeStateRec = outRec2;
			else if (OutRec1RightOfOutRec2(outRec2, outRec1)) holeStateRec = outRec1;
			else holeStateRec = GetLowermostRec(outRec1, outRec2);

			if (!JoinPoints(join, outRec1, outRec2)) continue;

			if (outRec1 == outRec2)
			{
				//instead of joining two polygons, we've just created a new one by
				//splitting one polygon into two.
				outRec1->Pts = join->OutPt1;
				outRec1->BottomPt = 0;
				outRec2 = CreateOutRec();
				outRec2->Pts = join->OutPt2;

				//update all OutRec2.Pts Idx's ...
				UpdateOutPtIdxs(*outRec2);

				if (Poly2ContainsPoly1(outRec2->Pts, outRec1->Pts))
				{
					//outRec1 contains outRec2 ...
					outRec2->IsHole = !outRec1->IsHole;
					outRec2->FirstLeft = outRec1;

					if (m_UsingPolyTree) FixupFirstLefts2(outRec2, outRec1);

					if ((outRec2->IsHole ^ m_ReverseOutput) == (Area(*outRec2) > 0))
						ReversePolyPtLinks(outRec2->Pts);

				}
				else if (Poly2ContainsPoly1(outRec1->Pts, outRec2->Pts))
				{
					//outRec2 contains outRec1 ...
					outRec2->IsHole = outRec1->IsHole;
					outRec1->IsHole = !outRec2->IsHole;
					outRec2->FirstLeft = outRec1->FirstLeft;
					outRec1->FirstLeft = outRec2;

					if (m_UsingPolyTree) FixupFirstLefts2(outRec1, outRec2);

					if ((outRec1->IsHole ^ m_ReverseOutput) == (Area(*outRec1) > 0))
						ReversePolyPtLinks(outRec1->Pts);
				}
				else
				{
					//the 2 polygons are completely separate ...
					outRec2->IsHole = outRec1->IsHole;
					outRec2->FirstLeft = outRec1->FirstLeft;

					//fixup FirstLeft pointers that may need reassigning to OutRec2
					if (m_UsingPolyTree) FixupFirstLefts1(outRec1, outRec2);
				}

			}
			else
			{
				//joined 2 polygons together ...

				outRec2->Pts = 0;
				outRec2->BottomPt = 0;
				outRec2->Idx = outRec1->Idx;

				outRec1->IsHole = holeStateRec->IsHole;
				if (holeStateRec == outRec2)
					outRec1->FirstLeft = outRec2->FirstLeft;
				outRec2->FirstLeft = outRec1;

				if (m_UsingPolyTree) FixupFirstLefts3(outRec2, outRec1);
			}
		}
	}

	ClipperOffsetEx::ClipperOffsetEx(double miterLimit, double arcTolerance)
	{
		this->MiterLimit = miterLimit;
		this->ArcTolerance = arcTolerance;
		m_lowest.X = -1;
	}
	//------------------------------------------------------------------------------

	ClipperOffsetEx::~ClipperOffsetEx()
	{
		Clear();
	}
	//------------------------------------------------------------------------------

	void ClipperOffsetEx::Clear()
	{
		for (int i = 0; i < m_polyNodes.ChildCount(); ++i)
			delete m_polyNodes.Childs[i];
		m_polyNodes.Childs.clear();
		m_lowest.X = -1;
	}
	//------------------------------------------------------------------------------

	void ClipperOffsetEx::AddPath(const Path& path, JoinType joinType, EndType endType)
	{
		int highI = (int)path.size() - 1;
		if (highI < 0) return;
		PolyNode* newNode = new PolyNode();
		newNode->m_jointype = joinType;
		newNode->m_endtype = endType;

		//strip duplicate points from path and also get index to the lowest point ...
		if (endType == etClosedLine || endType == etClosedPolygon)
			while (highI > 0 && path[0] == path[highI]) highI--;
		newNode->Contour.reserve(highI + 1);
		newNode->Contour.push_back(path[0]);
		int j = 0, k = 0;
		for (int i = 1; i <= highI; i++)
			if (newNode->Contour[j] != path[i])
			{
				j++;
				newNode->Contour.push_back(path[i]);
				if (path[i].Y > newNode->Contour[k].Y ||
					(path[i].Y == newNode->Contour[k].Y &&
						path[i].X < newNode->Contour[k].X)) k = j;
			}
		if (endType == etClosedPolygon && j < 2)
		{
			delete newNode;
			return;
		}
		m_polyNodes.AddChild(*newNode);

		//if this path's lowest pt is lower than all the others then update m_lowest
		if (endType != etClosedPolygon) return;
		if (m_lowest.X < 0)
			m_lowest = IntPoint(m_polyNodes.ChildCount() - 1, k);
		else
		{
			IntPoint ip = m_polyNodes.Childs[(int)m_lowest.X]->Contour[(int)m_lowest.Y];
			if (newNode->Contour[k].Y > ip.Y ||
				(newNode->Contour[k].Y == ip.Y &&
					newNode->Contour[k].X < ip.X))
				m_lowest = IntPoint(m_polyNodes.ChildCount() - 1, k);
		}
	}
	//------------------------------------------------------------------------------

	void ClipperOffsetEx::AddPaths(const Paths& paths, JoinType joinType, EndType endType)
	{
		for (Paths::size_type i = 0; i < paths.size(); ++i)
			AddPath(paths[i], joinType, endType);
	}
	//------------------------------------------------------------------------------

	void ClipperOffsetEx::FixOrientations()
	{
		//fixup orientations of all closed paths if the orientation of the
		//closed path with the lowermost vertex is wrong ...
		if (m_lowest.X >= 0 &&
			!Orientation(m_polyNodes.Childs[(int)m_lowest.X]->Contour))
		{
			for (int i = 0; i < m_polyNodes.ChildCount(); ++i)
			{
				PolyNode& node = *m_polyNodes.Childs[i];
				if (node.m_endtype == etClosedPolygon ||
					(node.m_endtype == etClosedLine && Orientation(node.Contour)))
					ReversePath(node.Contour);
			}
		}
		else
		{
			for (int i = 0; i < m_polyNodes.ChildCount(); ++i)
			{
				PolyNode& node = *m_polyNodes.Childs[i];
				if (node.m_endtype == etClosedLine && !Orientation(node.Contour))
					ReversePath(node.Contour);
			}
		}
	}
	//------------------------------------------------------------------------------

	void ClipperOffsetEx::ExecuteConst(Paths& solution, double delta, int step)
	{
		solution.clear();
		FixOrientations();
		DoConstOffset(delta, step);

		solution = m_destPolys;
	}

	void ClipperOffsetEx::DoConstOffset(double delta, int in_step)
	{
		m_destPolys.clear();
		m_delta = delta;

		//see offset_triginometry3.svg in the documentation folder ...
		if (MiterLimit > 2) m_miterLim = 2 / (MiterLimit * MiterLimit);
		else m_miterLim = 0.5;

		double y;
		if (ArcTolerance <= 0.0) y = def_arc_tolerance;
		else if (ArcTolerance > std::fabs(delta) * def_arc_tolerance)
			y = std::fabs(delta) * def_arc_tolerance;
		else y = ArcTolerance;
		//see offset_triginometry2.svg in the documentation folder ...
		double steps = pi / std::acos(1 - y / std::fabs(delta));
		if (steps > std::fabs(delta) * pi)
			steps = std::fabs(delta) * pi;  //ie excessive precision check

		if (in_step > 0) steps = (double)in_step;

		m_sin = std::sin(two_pi / steps);
		m_cos = std::cos(two_pi / steps);
		m_StepsPerRad = steps / two_pi;
		if (delta < 0.0) m_sin = -m_sin;

		m_destPolys.reserve(m_polyNodes.ChildCount() * 2);
		for (int i = 0; i < m_polyNodes.ChildCount(); i++)
		{
			PolyNode& node = *m_polyNodes.Childs[i];
			m_srcPoly = node.Contour;

			int len = (int)m_srcPoly.size();
			if (len == 0 || (delta <= 0 && (len < 3 || node.m_endtype != etClosedPolygon)))
				continue;

			m_destPoly.clear();
			if (len == 1)
			{
				if (node.m_jointype == jtRound)
				{
					double X = 1.0, Y = 0.0;
					for (cInt j = 1; j <= steps; j++)
					{
						m_destPoly.push_back(IntPoint(
							Round(m_srcPoly[0].X + X * delta),
							Round(m_srcPoly[0].Y + Y * delta)));
						double X2 = X;
						X = X * m_cos - m_sin * Y;
						Y = X2 * m_sin + Y * m_cos;
					}
				}
				else
				{
					double X = -1.0, Y = -1.0;
					for (int j = 0; j < 4; ++j)
					{
						m_destPoly.push_back(IntPoint(
							Round(m_srcPoly[0].X + X * delta),
							Round(m_srcPoly[0].Y + Y * delta)));
						if (X < 0) X = 1;
						else if (Y < 0) Y = 1;
						else X = -1;
					}
				}
				m_destPolys.push_back(m_destPoly);
				continue;
			}
			//build m_normals ...
			m_normals.clear();
			m_normals.reserve(len);
			for (int j = 0; j < len - 1; ++j)
				m_normals.push_back(GetUnitNormal(m_srcPoly[j], m_srcPoly[j + 1]));
			if (node.m_endtype == etClosedLine || node.m_endtype == etClosedPolygon)
				m_normals.push_back(GetUnitNormal(m_srcPoly[len - 1], m_srcPoly[0]));
			else
				m_normals.push_back(DoublePoint(m_normals[len - 2]));

			if (node.m_endtype == etClosedPolygon)
			{
				int k = len - 1;
				for (int j = 0; j < len; ++j)
					OffsetPoint(j, k, node.m_jointype);
				m_destPolys.push_back(m_destPoly);
			}
			else if (node.m_endtype == etClosedLine)
			{
				int k = len - 1;
				for (int j = 0; j < len; ++j)
					OffsetPoint(j, k, node.m_jointype);
				m_destPolys.push_back(m_destPoly);
				m_destPoly.clear();
				//re-build m_normals ...
				DoublePoint n = m_normals[len - 1];
				for (int j = len - 1; j > 0; j--)
					m_normals[j] = DoublePoint(-m_normals[j - 1].X, -m_normals[j - 1].Y);
				m_normals[0] = DoublePoint(-n.X, -n.Y);
				k = 0;
				for (int j = len - 1; j >= 0; j--)
					OffsetPoint(j, k, node.m_jointype);
				m_destPolys.push_back(m_destPoly);
			}
			else
			{
				int k = 0;
				for (int j = 1; j < len - 1; ++j)
					OffsetPoint(j, k, node.m_jointype);

				IntPoint pt1;
				if (node.m_endtype == etOpenButt)
				{
					int j = len - 1;
					pt1 = IntPoint((cInt)Round(m_srcPoly[j].X + m_normals[j].X *
						delta), (cInt)Round(m_srcPoly[j].Y + m_normals[j].Y * delta));
					m_destPoly.push_back(pt1);
					pt1 = IntPoint((cInt)Round(m_srcPoly[j].X - m_normals[j].X *
						delta), (cInt)Round(m_srcPoly[j].Y - m_normals[j].Y * delta));
					m_destPoly.push_back(pt1);
				}
				else
				{
					int j = len - 1;
					k = len - 2;
					m_sinA = 0;
					m_normals[j] = DoublePoint(-m_normals[j].X, -m_normals[j].Y);
					if (node.m_endtype == etOpenSquare)
						DoSquare(j, k);
					else
						DoRound(j, k);
				}

				//re-build m_normals ...
				for (int j = len - 1; j > 0; j--)
					m_normals[j] = DoublePoint(-m_normals[j - 1].X, -m_normals[j - 1].Y);
				m_normals[0] = DoublePoint(-m_normals[1].X, -m_normals[1].Y);

				k = len - 1;
				for (int j = k - 1; j > 0; --j) OffsetPoint(j, k, node.m_jointype);

				if (node.m_endtype == etOpenButt)
				{
					pt1 = IntPoint((cInt)Round(m_srcPoly[0].X - m_normals[0].X * delta),
						(cInt)Round(m_srcPoly[0].Y - m_normals[0].Y * delta));
					m_destPoly.push_back(pt1);
					pt1 = IntPoint((cInt)Round(m_srcPoly[0].X + m_normals[0].X * delta),
						(cInt)Round(m_srcPoly[0].Y + m_normals[0].Y * delta));
					m_destPoly.push_back(pt1);
				}
				else
				{
					k = 1;
					m_sinA = 0;
					if (node.m_endtype == etOpenSquare)
						DoSquare(0, 1);
					else
						DoRound(0, 1);
				}
				m_destPolys.push_back(m_destPoly);
			}
		}
	}

	void ClipperOffsetEx::Execute(Paths& solution, double delta)
	{
		solution.clear();
		FixOrientations();
		DoOffset(delta, -1);

		//now clean up 'corners' ...
		Clipper clpr;
		clpr.AddPaths(m_destPolys, ptSubject, true);
		if (delta > 0)
		{
			clpr.Execute(ctUnion, solution, pftPositive, pftPositive);
		}
		else
		{
			IntRect r = clpr.GetBounds();
			Path outer(4);
			outer[0] = IntPoint(r.left - 10, r.bottom + 10);
			outer[1] = IntPoint(r.right + 10, r.bottom + 10);
			outer[2] = IntPoint(r.right + 10, r.top - 10);
			outer[3] = IntPoint(r.left - 10, r.top - 10);

			clpr.AddPath(outer, ptSubject, true);
			clpr.ReverseSolution(true);
			clpr.Execute(ctUnion, solution, pftNegative, pftNegative);
			if (solution.size() > 0) solution.erase(solution.begin());
		}
	}
	//------------------------------------------------------------------------------

	void ClipperOffsetEx::Execute(PolyTree& solution, double delta)
	{
		solution.Clear();
		FixOrientations();
		DoOffset(delta, -1);

		//now clean up 'corners' ...
		Clipper clpr;
		clpr.AddPaths(m_destPolys, ptSubject, true);
		if (delta > 0)
		{
			clpr.Execute(ctUnion, solution, pftPositive, pftPositive);
		}
		else
		{
			IntRect r = clpr.GetBounds();
			Path outer(4);
			outer[0] = IntPoint(r.left - 10, r.bottom + 10);
			outer[1] = IntPoint(r.right + 10, r.bottom + 10);
			outer[2] = IntPoint(r.right + 10, r.top - 10);
			outer[3] = IntPoint(r.left - 10, r.top - 10);

			clpr.AddPath(outer, ptSubject, true);
			clpr.ReverseSolution(true);
			clpr.Execute(ctUnion, solution, pftNegative, pftNegative);
			//remove the outer PolyNode rectangle ...
			if (solution.ChildCount() == 1 && solution.Childs[0]->ChildCount() > 0)
			{
				PolyNode* outerNode = solution.Childs[0];
				solution.Childs.reserve(outerNode->ChildCount());
				solution.Childs[0] = outerNode->Childs[0];
				solution.Childs[0]->Parent = outerNode->Parent;
				for (int i = 1; i < outerNode->ChildCount(); ++i)
					solution.AddChild(*outerNode->Childs[i]);
			}
			else
				solution.Clear();
		}
	}
	//------------------------------------------------------------------------------

	void ClipperOffsetEx::DoOffset(double delta, int in_step)
	{
		m_destPolys.clear();
		m_delta = delta;

		//if Zero offset, just copy any CLOSED polygons to m_p and return ...
		if (NEAR_ZERO(delta))
		{
			m_destPolys.reserve(m_polyNodes.ChildCount());
			for (int i = 0; i < m_polyNodes.ChildCount(); i++)
			{
				PolyNode& node = *m_polyNodes.Childs[i];
				if (node.m_endtype == etClosedPolygon)
					m_destPolys.push_back(node.Contour);
			}
			return;
		}

		//see offset_triginometry3.svg in the documentation folder ...
		if (MiterLimit > 2) m_miterLim = 2 / (MiterLimit * MiterLimit);
		else m_miterLim = 0.5;

		double y;
		if (ArcTolerance <= 0.0) y = def_arc_tolerance;
		else if (ArcTolerance > std::fabs(delta) * def_arc_tolerance)
			y = std::fabs(delta) * def_arc_tolerance;
		else y = ArcTolerance;
		//see offset_triginometry2.svg in the documentation folder ...
		double steps = pi / std::acos(1 - y / std::fabs(delta));
		if (steps > std::fabs(delta) * pi)
			steps = std::fabs(delta) * pi;  //ie excessive precision check

		if (in_step > 0) steps = (double)in_step;

		m_sin = std::sin(two_pi / steps);
		m_cos = std::cos(two_pi / steps);
		m_StepsPerRad = steps / two_pi;
		if (delta < 0.0) m_sin = -m_sin;

		m_destPolys.reserve(m_polyNodes.ChildCount() * 2);
		for (int i = 0; i < m_polyNodes.ChildCount(); i++)
		{
			PolyNode& node = *m_polyNodes.Childs[i];
			m_srcPoly = node.Contour;

			int len = (int)m_srcPoly.size();
			if (len == 0 || (delta <= 0 && (len < 3 || node.m_endtype != etClosedPolygon)))
				continue;

			m_destPoly.clear();
			if (len == 1)
			{
				if (node.m_jointype == jtRound)
				{
					double X = 1.0, Y = 0.0;
					for (cInt j = 1; j <= steps; j++)
					{
						m_destPoly.push_back(IntPoint(
							Round(m_srcPoly[0].X + X * delta),
							Round(m_srcPoly[0].Y + Y * delta)));
						double X2 = X;
						X = X * m_cos - m_sin * Y;
						Y = X2 * m_sin + Y * m_cos;
					}
				}
				else
				{
					double X = -1.0, Y = -1.0;
					for (int j = 0; j < 4; ++j)
					{
						m_destPoly.push_back(IntPoint(
							Round(m_srcPoly[0].X + X * delta),
							Round(m_srcPoly[0].Y + Y * delta)));
						if (X < 0) X = 1;
						else if (Y < 0) Y = 1;
						else X = -1;
					}
				}
				m_destPolys.push_back(m_destPoly);
				continue;
			}
			//build m_normals ...
			m_normals.clear();
			m_normals.reserve(len);
			for (int j = 0; j < len - 1; ++j)
				m_normals.push_back(GetUnitNormal(m_srcPoly[j], m_srcPoly[j + 1]));
			if (node.m_endtype == etClosedLine || node.m_endtype == etClosedPolygon)
				m_normals.push_back(GetUnitNormal(m_srcPoly[len - 1], m_srcPoly[0]));
			else
				m_normals.push_back(DoublePoint(m_normals[len - 2]));

			if (node.m_endtype == etClosedPolygon)
			{
				int k = len - 1;
				for (int j = 0; j < len; ++j)
					OffsetPoint(j, k, node.m_jointype);
				m_destPolys.push_back(m_destPoly);
			}
			else if (node.m_endtype == etClosedLine)
			{
				int k = len - 1;
				for (int j = 0; j < len; ++j)
					OffsetPoint(j, k, node.m_jointype);
				m_destPolys.push_back(m_destPoly);
				m_destPoly.clear();
				//re-build m_normals ...
				DoublePoint n = m_normals[len - 1];
				for (int j = len - 1; j > 0; j--)
					m_normals[j] = DoublePoint(-m_normals[j - 1].X, -m_normals[j - 1].Y);
				m_normals[0] = DoublePoint(-n.X, -n.Y);
				k = 0;
				for (int j = len - 1; j >= 0; j--)
					OffsetPoint(j, k, node.m_jointype);
				m_destPolys.push_back(m_destPoly);
			}
			else
			{
				int k = 0;
				for (int j = 1; j < len - 1; ++j)
					OffsetPoint(j, k, node.m_jointype);

				IntPoint pt1;
				if (node.m_endtype == etOpenButt)
				{
					int j = len - 1;
					pt1 = IntPoint((cInt)Round(m_srcPoly[j].X + m_normals[j].X *
						delta), (cInt)Round(m_srcPoly[j].Y + m_normals[j].Y * delta));
					m_destPoly.push_back(pt1);
					pt1 = IntPoint((cInt)Round(m_srcPoly[j].X - m_normals[j].X *
						delta), (cInt)Round(m_srcPoly[j].Y - m_normals[j].Y * delta));
					m_destPoly.push_back(pt1);
				}
				else
				{
					int j = len - 1;
					k = len - 2;
					m_sinA = 0;
					m_normals[j] = DoublePoint(-m_normals[j].X, -m_normals[j].Y);
					if (node.m_endtype == etOpenSquare)
						DoSquare(j, k);
					else
						DoRound(j, k);
				}

				//re-build m_normals ...
				for (int j = len - 1; j > 0; j--)
					m_normals[j] = DoublePoint(-m_normals[j - 1].X, -m_normals[j - 1].Y);
				m_normals[0] = DoublePoint(-m_normals[1].X, -m_normals[1].Y);

				k = len - 1;
				for (int j = k - 1; j > 0; --j) OffsetPoint(j, k, node.m_jointype);

				if (node.m_endtype == etOpenButt)
				{
					pt1 = IntPoint((cInt)Round(m_srcPoly[0].X - m_normals[0].X * delta),
						(cInt)Round(m_srcPoly[0].Y - m_normals[0].Y * delta));
					m_destPoly.push_back(pt1);
					pt1 = IntPoint((cInt)Round(m_srcPoly[0].X + m_normals[0].X * delta),
						(cInt)Round(m_srcPoly[0].Y + m_normals[0].Y * delta));
					m_destPoly.push_back(pt1);
				}
				else
				{
					k = 1;
					m_sinA = 0;
					if (node.m_endtype == etOpenSquare)
						DoSquare(0, 1);
					else
						DoRound(0, 1);
				}
				m_destPolys.push_back(m_destPoly);
			}
		}
	}
	//------------------------------------------------------------------------------

	void ClipperOffsetEx::OffsetPoint(int j, int& k, JoinType jointype)
	{
		//cross product ...
		m_sinA = (m_normals[k].X * m_normals[j].Y - m_normals[j].X * m_normals[k].Y);
		if (std::fabs(m_sinA * m_delta) < 1.0)
		{
			//dot product ...
			double cosA = (m_normals[k].X * m_normals[j].X + m_normals[j].Y * m_normals[k].Y);
			if (cosA > 0) // angle => 0 degrees
			{
				m_destPoly.push_back(IntPoint(Round(m_srcPoly[j].X + m_normals[k].X * m_delta),
					Round(m_srcPoly[j].Y + m_normals[k].Y * m_delta)));
				return;
			}
			//else angle => 180 degrees   
		}
		else if (m_sinA > 1.0) m_sinA = 1.0;
		else if (m_sinA < -1.0) m_sinA = -1.0;

		if (m_sinA * m_delta < 0)
		{
			m_destPoly.push_back(IntPoint(Round(m_srcPoly[j].X + m_normals[k].X * m_delta),
				Round(m_srcPoly[j].Y + m_normals[k].Y * m_delta)));
			m_destPoly.push_back(m_srcPoly[j]);
			m_destPoly.push_back(IntPoint(Round(m_srcPoly[j].X + m_normals[j].X * m_delta),
				Round(m_srcPoly[j].Y + m_normals[j].Y * m_delta)));
		}
		else
			switch (jointype)
			{
			case jtMiter:
			{
				double r = 1 + (m_normals[j].X * m_normals[k].X +
					m_normals[j].Y * m_normals[k].Y);
				if (r >= m_miterLim) DoMiter(j, k, r); else DoSquare(j, k);
				break;
			}
			case jtSquare: DoSquare(j, k); break;
			case jtRound: DoRound(j, k); break;
			}
		k = j;
	}
	//------------------------------------------------------------------------------

	void ClipperOffsetEx::DoSquare(int j, int k)
	{
		double dx = std::tan(std::atan2(m_sinA,
			m_normals[k].X * m_normals[j].X + m_normals[k].Y * m_normals[j].Y) / 4);
		m_destPoly.push_back(IntPoint(
			Round(m_srcPoly[j].X + m_delta * (m_normals[k].X - m_normals[k].Y * dx)),
			Round(m_srcPoly[j].Y + m_delta * (m_normals[k].Y + m_normals[k].X * dx))));
		m_destPoly.push_back(IntPoint(
			Round(m_srcPoly[j].X + m_delta * (m_normals[j].X + m_normals[j].Y * dx)),
			Round(m_srcPoly[j].Y + m_delta * (m_normals[j].Y - m_normals[j].X * dx))));
	}
	//------------------------------------------------------------------------------

	void ClipperOffsetEx::DoMiter(int j, int k, double r)
	{
		double q = m_delta / r;
		m_destPoly.push_back(IntPoint(Round(m_srcPoly[j].X + (m_normals[k].X + m_normals[j].X) * q),
			Round(m_srcPoly[j].Y + (m_normals[k].Y + m_normals[j].Y) * q)));
	}
	//------------------------------------------------------------------------------

	void ClipperOffsetEx::DoRound(int j, int k)
	{
		double a = std::atan2(m_sinA,
			m_normals[k].X * m_normals[j].X + m_normals[k].Y * m_normals[j].Y);
		int steps = std::max((int)Round(m_StepsPerRad * std::fabs(a)), 1);

		double X = m_normals[k].X, Y = m_normals[k].Y, X2;
		for (int i = 0; i < steps; ++i)
		{
			m_destPoly.push_back(IntPoint(
				Round(m_srcPoly[j].X + X * m_delta),
				Round(m_srcPoly[j].Y + Y * m_delta)));
			X2 = X;
			X = X * m_cos - m_sin * Y;
			Y = X2 * m_sin + Y * m_cos;
		}
		m_destPoly.push_back(IntPoint(
			Round(m_srcPoly[j].X + m_normals[j].X * m_delta),
			Round(m_srcPoly[j].Y + m_normals[j].Y * m_delta)));
	}
}