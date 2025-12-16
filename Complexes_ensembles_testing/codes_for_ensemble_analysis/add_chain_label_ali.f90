program add_chain_label_ali
	implicit none
	
	character (len = 150) :: pdb_id_CA_P, pdb_id_CA_P_ali
	character (len = 1) :: chain_id_fill, chain_id_void, alt_loc
	character (len = 3) :: atom_type
	character (len = 4) :: atom_label
	character (len = 3) :: res_type
	character (len = 54) :: string
	integer :: num_lines, n, atom_num, res_num, i, count
	double precision :: x, y, z
	logical :: exist

	call getarg(1, pdb_id_CA_P)
	call getarg(2, pdb_id_CA_P_ali)

	open(11,file = trim(pdb_id_CA_P)//".pdb",status = 'old')
	num_lines = 0
	do
	    read(11,'(a54)',end=10) string
	    num_lines = num_lines + 1
	end do
10  close(11)
	
	
	! format:ATOM    110  CD1 PHE A  15      27.005  -7.442  11.099  1.00 19.10 
200 format(a4,i7,2x,a3,a1,a3,1x,a1,i4,4x,3f8.3)
	

	open(12,file = trim(pdb_id_CA_P)//".pdb",status = 'old')
	open(13,file = trim(pdb_id_CA_P_ali)//".pdb",status = 'old')


	inquire(file=trim(pdb_id_CA_P_ali)//"_labeled.pdb",exist=exist)
	if (exist) then
		open(14,file = trim(pdb_id_CA_P_ali)//"_labeled.pdb",status = 'old')
		close(14,status='delete')
	end if
	open(15,file = trim(pdb_id_CA_P_ali)//"_labeled.pdb",status = 'new')
	
	do i = 1, num_lines
		read(12,200) atom_label, atom_num, atom_type, alt_loc, res_type, chain_id_fill, res_num, x, y, z
		read(13,200) atom_label, atom_num, atom_type, alt_loc, res_type, chain_id_void, res_num, x, y, z
		write(15,200) atom_label, atom_num, atom_type, alt_loc, res_type, chain_id_fill, res_num, x, y, z
	end do

	close(12)
	close(13)
	close(15)

	open(16,file = trim(pdb_id_CA_P_ali)//".pdb",status = 'old')
	close(16,status='delete')

end program add_chain_label_ali