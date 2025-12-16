program run_pca_symm_red
	implicit none

	integer :: n_CA_P, n_exp_conf, i, j, k
	integer :: iostat, num_atom, num_res, count, info
	double precision :: x, y, z
	character (len = 1) :: alt_loc, chain
	character (len = 4) :: pdb_name, label_atom
	character (len = 100) :: chain_id
	character (len = 10) :: n_exp_conf_str, n_CA_P_str
	character (len = 3) :: res_type
	character (len = 3) :: atom_type
	double precision, allocatable :: coord_matrix(:,:), red_coord_matrix(:,:)
	double precision, allocatable :: cov_matrix(:,:), sacr_cov_matrix(:,:)
	double precision, allocatable :: evals(:), evals_set(:)
	double precision, allocatable :: evecs(:,:), work(:)
	double precision, allocatable :: sorted_evals(:), sorted_evecs(:,:)
	double precision, allocatable :: variance(:)
	logical :: exist
	double precision, allocatable :: exp_projections(:,:)
	double precision, allocatable :: avg_struct(:,:)
	double precision, allocatable :: red_conf1_coord_matrix(:,:)
	character (len = 1) :: ave_ref_flag
	double precision :: vl, vu
	integer, allocatable :: iwork(:), ifail(:)
	integer :: m, n


	external dsyevx


	call getarg(1, n_CA_P_str)
	call getarg(2, n_exp_conf_str)
	call getarg(3, ave_ref_flag)

	read (n_CA_P_str, *) n_CA_P
	read (n_exp_conf_str, *) n_exp_conf


	!!! READING THE PDB CODES IN ENSEMBLE.TXT AND SAVING THEIR 3D COORDINATES FOR PCA !!!

	open(11,file='ensemble.txt',status = 'old',iostat = iostat)
	if (iostat.ne.0) then
		print *, "I couldn't open the file ensemble.txt or the file does not exist!"
		stop
	end if

	allocate(coord_matrix(n_exp_conf,3*n_CA_P))
	count = 0
	do i = 1, n_exp_conf
100 format(a4,1x,a100)
200 format(a4,i7,2x,a3,a1,a3,1x,a1,i4,4x,3f8.3)
		read(11,100) pdb_name, chain_id
		open(12+count,file=pdb_name//"_"//trim(chain_id)//"_CA_P_ali_labeled.pdb",status = 'old')
		do j = 1, n_CA_P
			read(12+count,200) label_atom, num_atom, atom_type, alt_loc, res_type, chain, num_res, x, y, z
			coord_matrix(i,3*j-2) = x
			coord_matrix(i,3*j-1) = y
			coord_matrix(i,3*j) = z
		end do
		close(12+count)
		count = count + 1
	end do
	close(11)


	!!! PCA CALCULATION !!!

	allocate(red_coord_matrix(n_exp_conf,3*n_CA_P))
	allocate(avg_struct(1,3*n_CA_P))
	allocate(red_conf1_coord_matrix(n_exp_conf,3*n_CA_P))

	avg_struct = 0.0d0
	do i = 1, 3*n_CA_P
		do k = 1, n_exp_conf
			avg_struct(1,i) = avg_struct(1,i) +	coord_matrix(k,i)/dble(n_exp_conf)
		end do	
	end do

	if (ave_ref_flag.eq."R") then
		avg_struct(1,:) = coord_matrix(1,:)
		print *, "PCA performed with respect to the reference structure"
	elseif (ave_ref_flag.eq."A") then
		print *, "PCA performed with respect to the average structure"
	else
		print *, "You have to specify R for reference structure or A for average structure for PCA calculation!" 
		stop
	end if

	do i = 1, n_exp_conf
		do j = 1, 3*n_CA_P
			red_coord_matrix(i,j) = coord_matrix(i,j) - avg_struct(1,j)
			red_conf1_coord_matrix(i,j) = coord_matrix(i,j) - coord_matrix(1,j)
		end do
	end do

	allocate(cov_matrix(3*n_CA_P,3*n_CA_P))

	do i = 1, 3*n_CA_P
		do j = 1, 3*n_CA_P
			do k = 1, n_exp_conf
				cov_matrix(i,j) = cov_matrix(i,j) + (red_coord_matrix(k,i)*red_coord_matrix(k,j))/dble(n_exp_conf)
			end do
		end do
	end do

	allocate(sacr_cov_matrix(3*n_CA_P,3*n_CA_P))
	sacr_cov_matrix = cov_matrix

	if ((n_exp_conf - 1).ge.10) then
		m = n_exp_conf - 1
	else
		m = 10
	end if
	n = 3*n_CA_P

	allocate(evals(n))
	allocate(evecs(n,m))
	allocate(work(8*n))
	allocate(iwork(5*n))
	allocate(ifail(n))
	allocate(evals_set(m))

	call dsyevx('V','I','L',n,sacr_cov_matrix,n,vl,vu,n-m+1,n,0.00000001d0,m,evals,evecs,n,work,8*n,iwork,ifail,info)
	
	evals_set = evals(1:m)

	allocate(sorted_evals(m))
	allocate(sorted_evecs(n,m))

	call sort_eigs_desc(n,m,evals_set,evecs,sorted_evals,sorted_evecs)
	

	allocate(variance(m))
	do i = 1, m
		variance(i) = (sorted_evals(i)/sum(sorted_evals))
	end do

	inquire(file='pc_variances.txt',exist=exist)
	if (exist) then
		open(12+n_exp_conf,file='pc_variances.txt',status = 'old')
		close(12+n_exp_conf,status='delete')
	end if
300 format(f6.2)
	open(13+n_exp_conf,file='pc_variances.txt',status='new')
	do i = 1, m
		write(13+n_exp_conf,300) variance(i)*1.d2
	end do
	close(13+n_exp_conf)


	inquire(file="PCs.txt",exist=exist)
	if (exist) then
		open(14+n_exp_conf,file="PCs.txt",status = 'old')
		close(14+n_exp_conf,status='delete')
	end if
400 format(10f10.5)
	open(15+n_exp_conf,file="PCs.txt",status='new')
	do i = 1, 3*n_CA_P
		write(15+n_exp_conf,400) sorted_evecs(i,1:10)
	end do
	close(15+n_exp_conf)


	!!! PROJECT EXP CONFORMATIONS ON PCA SPACE !!!

	allocate(exp_projections(n_exp_conf,10))

	exp_projections = 0.0d0
	do j = 1, 10
		do k = 2, n_exp_conf
			do i = 1, 3*n_CA_P
				exp_projections(k,j) = exp_projections(k,j) + sorted_evecs(i,j)*red_conf1_coord_matrix(k,i)
			end do  
		end do
	end do

500 format(10f9.2)
	inquire(file='exp_ensemble_proj_PC.txt',exist=exist)
	if (exist) then
		open(16+n_exp_conf,file='exp_ensemble_proj_PC.txt',status = 'old')
		close(16+n_exp_conf,status='delete')
	end if
	open(17+n_exp_conf,file='exp_ensemble_proj_PC.txt',status='new')
	do k = 1, n_exp_conf
		write(17+n_exp_conf,500) exp_projections(k,1:10)
	end do
	close(17+n_exp_conf)



	!! SAVE COORDINATES AVERAGE STRUCTURE !!

	if (ave_ref_flag.eq."A") then
		inquire(file='avg_struct_pca.pdb',exist=exist)
		if (exist) then
			open(20,file='avg_struct_pca.pdb',status = 'old')
			close(20,status='delete')
		end if
		open(21,file='avg_struct_pca.pdb',status = 'new')
	else
		inquire(file='ref_struct_pca.pdb',exist=exist)
		if (exist) then
			open(20,file='ref_struct_pca.pdb',status = 'old')
			close(20,status='delete')
		end if
		open(21,file='ref_struct_pca.pdb',status = 'new')
	end if

	open(11,file='ensemble.txt',status = 'old')
	read(11,100) pdb_name, chain_id
	close(11)

	open(22,file=pdb_name//"_"//trim(chain_id)//"_CA_P_ali_labeled.pdb",status = 'old')
	
	do i = 1, n_CA_P
		read(22,200) label_atom, num_atom, atom_type, alt_loc, res_type, chain, num_res, x, y, z
		x = avg_struct(1,3*i-2)
		y = avg_struct(1,3*i-1)
		z = avg_struct(1,3*i)
		write(21,200) label_atom, num_atom, atom_type, alt_loc, res_type, chain, num_res, x, y, z
	end do
	
	close(21)
	close(22)

end program run_pca_symm_red


subroutine sort_eigs_desc(n,m,evals,evecs,sorted_evals,sorted_evecs) 
	implicit none
	
	integer :: i, j
	integer, intent(in) :: n,m
	integer :: sort_index(m)
	double precision, intent(in) :: evals(m), evecs(n,m)
	double precision, intent(out) :: sorted_evals(n), sorted_evecs(n,m)
	double precision :: max_value, sacr_evals(m)


	!Define sort_indexes vector - evals in descending order!

	sacr_evals = evals
	do i = 1,m
		max_value = minval(sacr_evals)
		do j = 1,m
			if (sacr_evals(j) >= max_value) then
				max_value = sacr_evals(j)
				sort_index(i) = j
			end if
		end do
		sacr_evals(sort_index(i)) = minval(sacr_evals) - 1D+0
	end do

	!Sort evals and evecs!

	do i = 1,m
		sorted_evals(i) = evals(sort_index(i))
		do j = 1,n
			sorted_evecs(j,i) = evecs(j,sort_index(i))
		end do
	end do

end subroutine sort_eigs_desc