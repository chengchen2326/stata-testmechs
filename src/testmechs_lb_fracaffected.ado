program define testmechs_lb_fracaffected, rclass
    version 16.0

    syntax varlist(min=3 numeric) [if] [in] [, atgroup(string) numybins(string) maxdefiersshare(string) allowmindefiers]

    if ("`numybins'" == "") local numybins 5
    if ("`maxdefiersshare'" == "") local maxdefiersshare 0
    local numybins = real("`numybins'")
    local maxdefiersshare = real("`maxdefiersshare'")

    marksample touse
    local nvars : word count `varlist'
    local d : word 1 of `varlist'
    local y : word `nvars' of `varlist'
    local mvars
    forvalues i = 2/`=`nvars'-1' {
        local mi : word `i' of `varlist'
        local mvars `mvars' `mi'
    }

    if `maxdefiersshare' < 0 {
        di as err "maxdefiersshare() must be nonnegative"
        exit 198
    }

    tempvar wt touse2 ywork missm mgroup
    quietly gen double `wt' = 1 if `touse'
    quietly egen byte `missm' = rowmiss(`mvars') if `touse'
    quietly gen byte `touse2' = `touse' & !missing(`d', `y', `wt') & `missm' == 0
    quietly egen long `mgroup' = group(`mvars') if `touse2'

    quietly count if `touse2' & !inlist(`d', 0, 1)
    if r(N) > 0 {
        di as err "treatment variable must be coded 0/1"
        exit 198
    }

    quietly count if `touse2' & `d' == 1
    local n1 = r(N)
    quietly count if `touse2' & `d' == 0
    local n0 = r(N)
    if (`n1' == 0 | `n0' == 0) {
        di as err "both treatment groups (d==0 and d==1) must be present"
        exit 2000
    }

    if `numybins' > 0 quietly xtile `ywork' = `y' if `touse2', nq(`numybins')
    else quietly gen double `ywork' = `y' if `touse2'

    quietly levelsof `mgroup' if `touse2', local(mlevels)
    local K : word count `mlevels'
    if `K' < 2 {
        di as err "mediator must have at least 2 observed levels"
        exit 198
    }

    quietly summarize `wt' if `touse2' & `d' == 1, meanonly
    scalar W1 = r(sum)
    quietly summarize `wt' if `touse2' & `d' == 0, meanonly
    scalar W0 = r(sum)

    tempname p1 p0 maxdiff
    matrix `p1' = J(`K',1,0)
    matrix `p0' = J(`K',1,0)
    matrix `maxdiff' = J(`K',1,0)

    local k = 0
    foreach mv of local mlevels {
        local ++k
        quietly summarize `wt' if `touse2' & `d' == 1 & `mgroup' == `mv', meanonly
        matrix `p1'[`k',1] = r(sum) / W1
        quietly summarize `wt' if `touse2' & `d' == 0 & `mgroup' == `mv', meanonly
        matrix `p0'[`k',1] = r(sum) / W0

        quietly levelsof `ywork' if `touse2' & `mgroup' == `mv', local(yvals)
        scalar thisdiff = 0
        foreach yv of local yvals {
            quietly summarize `wt' if `touse2' & `d' == 1 & `mgroup' == `mv' & `ywork' == `yv', meanonly
            scalar p1y = r(sum) / W1
            quietly summarize `wt' if `touse2' & `d' == 0 & `mgroup' == `mv' & `ywork' == `yv', meanonly
            scalar p0y = r(sum) / W0
            scalar thisdiff = thisdiff + max(p1y - p0y, 0)
        }
        matrix `maxdiff'[`k',1] = thisdiff
    }

    local atindex 0
    if "`atgroup'" != "" {
        local atgroup_num = real("`atgroup'")
        local k = 0
        foreach mv of local mlevels {
            local ++k
            if (`mv' == `atgroup_num') local atindex `k'
        }
        if `atindex' == 0 {
            di as err "atgroup() must equal one observed mediator level"
            exit 198
        }
    }

    tempfile lpinput lpout pysrc
    file open fh using `lpinput', write text replace
    file write fh "K `K'" _n
    file write fh "atindex `atindex'" _n
    file write fh "maxdef `maxdefiersshare'" _n
    local allowmin = cond("`allowmindefiers'" != "", 1, 0)
    file write fh "allowmin `allowmin'" _n
    forvalues i = 1/`K' {
        local p1i = `p1'[`i',1]
        local p0i = `p0'[`i',1]
        local mdi = `maxdiff'[`i',1]
        file write fh "row `i' `p1i' `p0i' `mdi'" _n
    }
    file close fh

    file open py using `pysrc', write text replace
    file write py "import sys" _n
    file write py "EPS=1e-9" _n
    file write py "INF=1e100" _n
    file write py "class LPSolver:" _n
    file write py "    def __init__(self,A,b,c):" _n
    file write py "        self.m=len(b); self.n=len(c)" _n
    file write py "        self.N=list(range(self.n))+[-1]" _n
    file write py "        self.B=[self.n+i for i in range(self.m)]" _n
    file write py "        self.D=[[0.0]*(self.n+2) for _ in range(self.m+2)]" _n
    file write py "        for i in range(self.m):" _n
    file write py "            for j in range(self.n): self.D[i][j]=A[i][j]" _n
    file write py "            self.D[i][self.n]=-1.0" _n
    file write py "            self.D[i][self.n+1]=b[i]" _n
    file write py "        for j in range(self.n): self.D[self.m][j]=-c[j]" _n
    file write py "        self.D[self.m+1][self.n]=1.0" _n
    file write py "    def pivot(self,r,s):" _n
    file write py "        D=self.D; m=self.m; n=self.n" _n
    file write py "        inv=1.0/D[r][s]" _n
    file write py "        for i in range(m+2):" _n
    file write py "            if i==r: continue" _n
    file write py "            for j in range(n+2):" _n
    file write py "                if j==s: continue" _n
    file write py "                D[i][j]-=D[r][j]*D[i][s]*inv" _n
    file write py "        for j in range(n+2):" _n
    file write py "            if j!=s: D[r][j]*=inv" _n
    file write py "        for i in range(m+2):" _n
    file write py "            if i!=r: D[i][s]*=-inv" _n
    file write py "        D[r][s]=inv" _n
    file write py "        self.B[r],self.N[s]=self.N[s],self.B[r]" _n
    file write py "    def simplex(self,phase):" _n
    file write py "        x=self.m+1 if phase==1 else self.m" _n
    file write py "        while True:" _n
    file write py "            s=-1" _n
    file write py "            for j in range(self.n+1):" _n
    file write py "                if phase==2 and self.N[j]==-1: continue" _n
    file write py "                if s==-1 or self.D[x][j] < self.D[x][s]-EPS or (abs(self.D[x][j]-self.D[x][s])<=EPS and self.N[j]<self.N[s]): s=j" _n
    file write py "            if self.D[x][s] >= -EPS: return True" _n
    file write py "            r=-1" _n
    file write py "            for i in range(self.m):" _n
    file write py "                if self.D[i][s] <= EPS: continue" _n
    file write py "                if r==-1: r=i" _n
    file write py "                else:" _n
    file write py "                    lhs=self.D[i][self.n+1]/self.D[i][s]" _n
    file write py "                    rhs=self.D[r][self.n+1]/self.D[r][s]" _n
    file write py "                    if lhs < rhs-EPS or (abs(lhs-rhs)<=EPS and self.B[i] < self.B[r]): r=i" _n
    file write py "            if r==-1: return False" _n
    file write py "            self.pivot(r,s)" _n
    file write py "    def solve_max(self):" _n
    file write py "        r=0" _n
    file write py "        for i in range(1,self.m):" _n
    file write py "            if self.D[i][self.n+1] < self.D[r][self.n+1]: r=i" _n
    file write py "        if self.D[r][self.n+1] < -EPS:" _n
    file write py "            self.pivot(r,self.n)" _n
    file write py "            if (not self.simplex(1)) or self.D[self.m+1][self.n+1] < -EPS: return None,None,False" _n
    file write py "            if abs(self.D[self.m+1][self.n+1]) > EPS: return None,None,False" _n
    file write py "            rr=-1" _n
    file write py "            for i in range(self.m):" _n
    file write py "                if self.B[i]==-1: rr=i; break" _n
    file write py "            if rr!=-1:" _n
    file write py "                s=0" _n
    file write py "                for j in range(1,self.n+1):" _n
    file write py "                    if self.D[rr][j] < self.D[rr][s]-EPS or (abs(self.D[rr][j]-self.D[rr][s])<=EPS and self.N[j] < self.N[s]): s=j" _n
    file write py "                self.pivot(rr,s)" _n
    file write py "        if not self.simplex(2): return None,None,True" _n
    file write py "        x=[0.0]*self.n" _n
    file write py "        for i in range(self.m):" _n
    file write py "            if self.B[i] < self.n: x[self.B[i]] = self.D[i][self.n+1]" _n
    file write py "        return x,self.D[self.m][self.n+1],False" _n
    file write py "def solve_min(A,b,c):" _n
    file write py "    x,val,unb=LPSolver(A,b,[-v for v in c]).solve_max()" _n
    file write py "    if x is None: return None,None,unb" _n
    file write py "    return x,-val,unb" _n
    file write py "inp,outp=sys.argv[1],sys.argv[2]" _n
    file write py "lines=[x.strip().split() for x in open(inp) if x.strip()]" _n
    file write py "K=int(lines[0][1]); atindex=int(lines[1][1]); maxdef=float(lines[2][1]); allowmin=int(lines[3][1])" _n
    file write py "rows=lines[4:4+K]" _n
    file write py "p1=[float(r[2]) for r in rows]; p0=[float(r[3]) for r in rows]; md=[float(r[4]) for r in rows]" _n
    file write py "def idx(i,j): return i*K+j" _n
    file write py "n=K*K" _n
    file write py "A=[]; b=[]" _n
    file write py "for j in range(K):" _n
    file write py "    row=[0.0]*n" _n
    file write py "    for i in range(K): row[idx(i,j)] = 1.0" _n
    file write py "    A.append(row); b.append(p1[j]); A.append([-v for v in row]); b.append(-p1[j])" _n
    file write py "for i in range(K):" _n
    file write py "    row=[0.0]*n" _n
    file write py "    for j in range(K): row[idx(i,j)] = 1.0" _n
    file write py "    A.append(row); b.append(p0[i]); A.append([-v for v in row]); b.append(-p0[i])" _n
    file write py "c=[1.0 if i>j else 0.0 for i in range(K) for j in range(K)]" _n
    file write py "x,val,unb=solve_min(A,b,c)" _n
    file write py "if x is None: raise RuntimeError('feasibility LP failed')" _n
    file write py "min_def=val" _n
    file write py "if min_def>maxdef:" _n
    file write py "    if allowmin==1: maxdef=min_def+1e-6" _n
    file write py "    else: raise RuntimeError('data incompatible with maxdefiersshare when allowmindefiers is off')" _n
    file write py "groups=list(range(K)) if atindex==0 else [atindex-1]" _n
    file write py "N=n+K+1; tix=N-1" _n
    file write py "A=[]; b=[]" _n
    file write py "for j in range(K):" _n
    file write py "    row=[0.0]*N" _n
    file write py "    for i in range(K): row[idx(i,j)] = 1.0" _n
    file write py "    row[tix] = -p1[j]" _n
    file write py "    A.append(row); b.append(0.0); A.append([-v for v in row]); b.append(0.0)" _n
    file write py "for i in range(K):" _n
    file write py "    row=[0.0]*N" _n
    file write py "    for j in range(K): row[idx(i,j)] = 1.0" _n
    file write py "    row[tix] = -p0[i]" _n
    file write py "    A.append(row); b.append(0.0); A.append([-v for v in row]); b.append(0.0)" _n
    file write py "row=[0.0]*N" _n
    file write py "for g in groups: row[idx(g,g)] = 1.0" _n
    file write py "A.append(row); b.append(1.0); A.append([-v for v in row]); b.append(-1.0)" _n
    file write py "for k in range(K):" _n
    file write py "    row=[0.0]*N" _n
    file write py "    for i in range(K):" _n
    file write py "        if i!=k: row[idx(i,k)] = -1.0" _n
    file write py "    row[n+k] = -1.0" _n
    file write py "    row[tix] = md[k]" _n
    file write py "    A.append(row); b.append(0.0)" _n
    file write py "row=[0.0]*N" _n
    file write py "for i in range(K):" _n
    file write py "    for j in range(K):" _n
    file write py "        if i>j: row[idx(i,j)] = 1.0" _n
    file write py "row[tix] = -maxdef" _n
    file write py "A.append(row); b.append(0.0)" _n
    file write py "obj=[0.0]*N" _n
    file write py "for g in groups: obj[n+g] = 1.0" _n
    file write py "x2,val2,unb2=solve_min(A,b,obj)" _n
    file write py "if x2 is None: raise RuntimeError('fractional LP failed')" _n
    file write py "with open(outp,'w') as f:" _n
    file write py "    f.write('lb %s\n' % val2)" _n
    file write py "    f.write('min_defier_share %s\n' % min_def)" _n
    file write py "    f.write('maxdefiersshare_used %s\n' % maxdef)" _n
    file close py

    quietly shell python3 `pysrc' `lpinput' `lpout'
    if _rc {
        di as err "python-based LP solver failed"
        exit 498
    }

    tempname lb min_def maxdef_used
    scalar `lb' = .
    scalar `min_def' = .
    scalar `maxdef_used' = .

    file open rf using `lpout', read text
    file read rf line
    while r(eof)==0 {
        local key : word 1 of `line'
        local val : word 2 of `line'
        if "`key'" == "lb" scalar `lb' = real("`val'")
        if "`key'" == "min_defier_share" scalar `min_def' = real("`val'")
        if "`key'" == "maxdefiersshare_used" scalar `maxdef_used' = real("`val'")
        file read rf line
    }
    file close rf

    return scalar lb = `lb'
    return scalar min_defier_share = `min_def'
    return scalar maxdefiersshare_used = `maxdef_used'

    di as txt "testmechs_lb_fracaffected"
    di as res "  lower bound = " %9.6f `lb'
end
