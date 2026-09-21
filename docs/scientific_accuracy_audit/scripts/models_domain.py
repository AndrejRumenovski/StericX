"""Independent known-geometry applicability/interval and numerical edge checks."""
from models_audit import *
from scipy.spatial.distance import cdist

def main():
 requests=[];refs={};mean=np.array([10.,-5.]);scale=np.array([2.,20.]);points=np.array([[-2,0],[-1,-1],[-1,1],[0,-2],[0,0],[0,2],[1,-1],[1,1],[2,0]],float)
 def geometry(points,regularize=0):
  n=len(points);D=np.c_[np.ones(n),points];N=D.T@D;N[1:,1:]+=regularize*np.eye(points.shape[1]);inv=np.linalg.inv(N);dist=cdist(points,points);np.fill_diagonal(dist,np.inf);nn=dist.min(1)
  return {'feature_indices':[1,2],'means':mean.tolist(),'scales':scale.tolist(),'xtx_inverse':inv.tolist(),'observations':n,'parameters':3,'residual_standard_error':.5,'warning_leverage':9/n,'standardized_training_points':points.tolist(),'training_labels':[str(i) for i in range(n)],'neighbor_calibration':{'mean':float(nn.mean()),'standard_deviation':float(nn.std()),'median':float(np.median(nn)),'maximum':float(nn.max()),'threshold':float(nn.max()),'rule':'audit independent max nearest other'}}
 g=geometry(points);ranges=[{'feature':'L_boltz','minimum':6.,'maximum':14.},{'feature':'B1_boltz','minimum':-45.,'maximum':35.}]
 for name,z in [('training',points[1]),('interpolation',np.array([.2,.3])),('inside_box_outside_convex_hull',np.array([1.9,1.9])),('extreme',np.array([100.,-100.]))]:
  vals=np.r_[1,mean+scale*z,0,0,0,0,0];v={'op':'model','id':name,'geometry':g,'features':vals.tolist(),'ranges':ranges,'prediction':2.};requests.append(v)
  h=float(np.r_[1,z]@np.asarray(g['xtx_inverse'])@np.r_[1,z]);mah=float(np.sqrt(z@np.linalg.inv(np.cov(points,rowvar=False))@z));t=float(stats.t.ppf(.975,6));refs[name]={'leverage':h,'mahalanobis':mah,'nearest_distance':float(cdist(z[None],points).min()),'pi':[2-t*.5*np.sqrt(1+h),2+t*.5*np.sqrt(1+h)],'ci':[2-t*.5*np.sqrt(h),2+t*.5*np.sqrt(h)]}
 # A rank-deficient covariance, stabilized exactly as the fit's documented floor.
 singular=np.array([[t,t] for t in [-2,-1,0,0,1,2]],float);sg=geometry(singular,1e-10)
 requests.append({'op':'model','id':'singular_covariance_off_manifold','geometry':sg,'features':[1,12,-25,0,0,0,0,0],'ranges':ranges,'prediction':0});refs['singular_covariance_off_manifold']={'matrix_rank':int(np.linalg.matrix_rank(np.cov(singular,rowvar=False))),'expected':'Ordinary Mahalanobis covariance inverse is undefined. Ridge-stabilized surrogate must be labeled as such.'}
 for df in [1.,2.,3.,8.,29.,100.,1e6]:
  for alpha in [.05,.01,.001,1e-8]:
   name=f't_df{df}_a{alpha}';requests.append({'op':'model','id':name,'df':df,'alpha':alpha});refs[name]={'t_quantile':float(stats.t.isf(alpha/2,df))}
 inp=M/'inputs/domain_requests.jsonl';freeze(inp,(''.join(json.dumps(v)+'\n' for v in requests)).encode());out=M/'raw/domain_observer.jsonl'
 if not out.exists():
  p=subprocess.run([str(A/'frozen/bin/stericx-audit-observer')],input=inp.read_bytes(),capture_output=True);freeze(out,p.stdout);freeze(M/'raw/domain_observer.stderr',p.stderr)
 rows=[]
 for line in out.read_text().splitlines():
  v=json.loads(line);ref=refs[v['id']];row={'id':v['id'],'stericx':v,'independent':ref}
  if 't_quantile' in ref:row['absolute_error']=abs(v['t_quantile']-ref['t_quantile']);row['relative_error']=row['absolute_error']/ref['t_quantile']
  rows.append(row)
 save(M/'results/domain_comparison.json',rows)

if __name__=='__main__':main()
