import glob
import os
import sys
import get_atten_primary
import get_atten_scatter
import shutil

def reconstruct(primaries_per_proj,num_jobs):
    
    # compute attenuation with scatter and reconstruct
    proj_path_withscatter='./withscatter'
    get_atten_scatter.attenuation(primaries_per_proj,proj_path_withscatter,num_jobs)
    reconstruct_steps(proj_path_withscatter)

    # compute attenuation without scatter and reconstruct
    proj_path_withoutscatter='./withoutscatter'
    get_atten_primary.attenuation(proj_path_withoutscatter)
    reconstruct_steps(proj_path_withoutscatter)


    # copy real cbct, and reconstructed with and without scatter, and mask to cbct_images folder
    copy_directory = f'../../cbct_images/{os.path.basename(os.getcwd())}'

    if not os.path.exists(copy_directory):
        os.makedirs(copy_directory)

    filename = f'./withscatter/Reconstruction/fdk_rotated.mha'
    if os.path.exists(filename):
        shutil.copy(filename,f'{copy_directory}/fdk_rotated_withscatter.mha')

    filename = f'./withoutscatter/Reconstruction/fdk_rotated.mha'
    if os.path.exists(filename):
        shutil.copy(filename,f'{copy_directory}/fdk_rotated_withoutscatter.mha')

    filename = f'./withscatter/Reconstruction/mask.mha'
    if os.path.exists(filename):
        shutil.copy(filename,f'{copy_directory}/mask_withscatter.mha')

    #filename = f'./data/cbct.0.nii'
    #if os.path.exists(filename):
        #shutil.copy(filename,copy_directory)

    # get patient mask
    extractpatient_command = f'clitkExtractPatient -i ./data/CT.nii -o {copy_directory}/patient_mask.mha'
    os.system(extractpatient_command)


def reconstruct_steps(projections_path):
    reconstructed_img_path= f'{projections_path}/Reconstruction'
    if not os.path.exists(reconstructed_img_path):
        os.makedirs(reconstructed_img_path)

    # reconstruct image 
    fdkcommand = f'rtkfdk -p {projections_path} -r .*.mha --pad=0.1 -o {reconstructed_img_path}/fdk.mha -g data/elektaGeometry.xml --dimension 410,264,410'
    os.system(fdkcommand)
    fovcommand = f'rtkfieldofview --geometry data/elektaGeometry.xml -d -p {projections_path} -r .*.mha --reconstruction {reconstructed_img_path}/fdk.mha  --output {reconstructed_img_path}/fdk.mha --verbose'
    #--hardware=cuda
    os.system(fovcommand)
    # save the mask
    maskcommand = f'rtkfieldofview --geometry data/elektaGeometry.xml -d -p {projections_path} -r .*.mha --reconstruction {reconstructed_img_path}/fdk.mha  --output {reconstructed_img_path}/mask.mha -m'
    os.system(maskcommand)

    # rotate reconstructed image and mask
    rotatecommand = f'clitkAffineTransform -i {reconstructed_img_path}/fdk.mha -o {reconstructed_img_path}/fdk_rotated.mha -m data/matriceRTK2CBCT.mat --pad=0 --transform_grid'
    os.system(rotatecommand)
    rotatemaskcommand = f'clitkAffineTransform -i {reconstructed_img_path}/mask.mha -o {reconstructed_img_path}/mask.mha -m data/matriceRTK2CBCT.mat --pad=0 --transform_grid'
    os.system(rotatemaskcommand)

    # rescale reconstructed image
    #rescale_pixels.rescale(f'{reconstructed_img_path}/fdk_rotated.mha')


if __name__ == "__main__":
    prim_per_proj = int(sys.argv[1])
    num_jobs = int(sys.argv[2])
    reconstruct(prim_per_proj,num_jobs)

