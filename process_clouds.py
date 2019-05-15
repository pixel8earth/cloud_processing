import sys, os, json 
import pdal
import numpy as np
from pathlib import Path

try:
    import trimesh
    from trimesh.registration import icp
    has_trimesh = True
except:
    has_trimesh = False

def execute(pipeline):
  pipeline = pdal.Pipeline(json.dumps(pipeline))
  pipeline.validate()
  pipeline.execute()

def reproject(model, output, proj="EPSG:32613"):
  pipeline = {
    "pipeline": [
        {
          "type" : "readers.text",
          "filename" : model
        }, {
          "type": "filters.reprojection",
          "in_srs": "+proj=geocent +ellps=WGS84 +datum=WGS84",
          "out_srs": proj
        },{
          "type" : "writers.text",
          "order":"X:5, Y:5, Z:5, Red:0, Green:0, Blue:0", 
          "keep_unspecified":"false",
          "quote_header": "false",
          "filename" : output
        }
    ]
  }
  execute(pipeline)

def newply(name, vals, header):
    if len(vals):
        with open(name, "w") as f:
            header[2] = 'element vertex {}\n'.format(len(vals))
            [f.write(h) for h in header]
            for v in vals:
                f.write(v)

def split_cams(model_path):
    """ 
      Splits the camera x,y,z lines from ply files 
      overwrites the given model path and creates the given path
    """
    cam_path = model_path.replace('.ply', '-cams.ply') 
    if os.path.exists(model_path):
        cams, pnts, header = [], [], []
        with open(model_path, "r") as f:
            for l in f.readlines():
                try:
                    vals = list(map(float, l.strip("\n").split(' ')))
                    x,y,z,r,g,b = vals 
                    if [r,g,b] == [0.0, 255.0, 0.0]:
                        cams.append(l)
                    else:
                        pnts.append(l)
                except Exception as err: 
                    header.append(l)
        newply(model_path, pnts, header)
        newply(cam_path, cams, header)

def apply_transform(csv, tfm, output):
    """ Applies a 4x4 transformation to a point cloud """ 
    tfm_str = ' '.join(list(map(str, tfm.flatten())))
    pipeline = [
      { "type": "readers.text", "filename": csv},
      {
        "type":"filters.transformation",
        "matrix": tfm_str
      },
      { 
        "type":"writers.text", 
        "order":"X:5, Y:5, Z:5, Red:0, Green:0, Blue:0",
        "keep_unspecified":"false",
        "quote_header": "false",
        "filename": output 
      }
    ]
    execute(pipeline)

def remove_outliers_pdal(model, k=8, mult=0.5, **kwargs):
    output = model.replace('.ply', '-clean.csv')
    pipeline = {
      "pipeline": [
          {
            "type" : "readers.ply",
            "filename" : model
          },{
            "type":"filters.elm"
          },{
            "type":"filters.outlier",
            "method": "statistical",
            "mean_k": k,
            "multiplier": mult
          },{
            "type": "filters.range",
            "limits": "Classification![7:7]",
          },{
            "type" : "writers.text",
            "order":"X:5, Y:5, Z:5, Red:0, Green:0, Blue:0",
            "keep_unspecified":"false",
            "quote_header": "false",
            "filename" : output
          }
        ]
    }
    execute(pipeline)

def clean(model, **kwargs):
    split_cams(str(model))
    remove_outliers_pdal(str(model), **kwargs)

def perform_icp(projected, output):
    msh = trimesh.load_mesh("poisson-utm.ply")
    with open(projected) as f:
        pts = np.asarray([[float(e) for e in row.strip().split(",")] for row in f.readlines()[1:]])[:,:3]
    pts[:,2] = pts[:,2] + np.mean(msh.vertices[:,2]) - np.mean(pts[:,2])
    z_extents = msh.bounds[:,2]
    z_levels  = np.arange(*z_extents, step=0.25)
    sections = msh.section_multiplane(plane_origin=[0,0,0],
                                       plane_normal=[0,0,1],
                                       heights=z_levels)
    ref_pts = np.vstack([np.hstack([section.vertices, z*np.ones((section.vertices.shape[0], 1))])
                                       for section, z in zip(sections[1:], z_levels[1:])])

    tfm, pts_tfm, _ = icp(pts, ref_pts, reflection=False, translation=True, scale=False)
    with open(output.replace('.csv', '-trimesh.csv'), 'w') as fh:
        for p in pts_tfm:
            fh.write(','.join(list(map(str,p))) + '\n')
            
    apply_transform(projected, tfm, transformed)
    return tfm

def perform_icp2(projected, output):
    msh = trimesh.load_mesh("poisson-utm.ply")
    with open(projected) as f:
        pts = np.asarray([[float(e) for e in row.strip().split(",")] for row in f.readlines()[1:]])[:,:3]
    tfm, pts_tfm, _ = icp(pts, msh, reflection=False, translation=True, scale=False)
    with open(output.replace('.csv', '-trimesh.csv'), 'w') as fh:
        for p in pts_tfm:
            fh.write(','.join(list(map(str,p))) + '\n')

    apply_transform(projected, tfm, transformed)
    return tfm


if __name__ == "__main__":
    stream = sys.argv[1]

    root = Path(f'./{stream}')
    model = str(root / 'model.ply')
    geo_model = str(root / 'geo-model.ply')
    geo_clean = str(root / 'geo-model-clean.csv')
    projected = str(root / 'geo-model-utm.csv')
    geo_transform = str(root / 'geo-transform.npy')
    transformed = str(root / 'geo-model-transformed.csv')
    temp = str(root / 'sfm.csv')

    if os.path.exists(model):
        clean(model)

    if os.path.exists(geo_model):
        clean(geo_model, k=8, mult=1.5)
        # reproject from ecef to UTM
        reproject(geo_clean, projected)

        # Mesh based ICP
        if has_trimesh:
            tfm = perform_icp(temp, transformed)
            #tfm = perform_icp2(projected, transformed)
            np.save(geo_transform, tfm)
        else:
            print('No trimesh... what to do')
    print('Finished Post Processing')
