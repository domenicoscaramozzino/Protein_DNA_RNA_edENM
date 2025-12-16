program write_CA_P_C2_C1
	implicit none
	
	character (len = 20) :: pdb_id
	character (len = 100) :: chain_to_read
	character (len = 1) :: chain_id, alt_loc
	character (len = 3) :: atom_type
	character (len = 4) :: atom_label
	character (len = 3) :: res_type
	character (len = 80) :: string
	integer :: pdb_iostat, num_atoms, n, atom_num, res_num, num_chains, i, count
	double precision :: x, y, z
	logical :: exist
	
	call getarg(1, pdb_id)
	call getarg(2, chain_to_read)

	open(11,file = trim(pdb_id)//".pdb",status = 'old',iostat = pdb_iostat)
	if (pdb_iostat.ne.0) then
		print *, "I couldn't open the file ", trim(pdb_id)//".pdb" ," or the file does not exist!"
		stop
	end if

	inquire(file=trim(pdb_id)//"_ATOM.pdb",exist=exist)
	if (exist) then
		open(12,file = trim(pdb_id)//"_ATOM.pdb",status = 'old')
		close(12,status='delete')
	end if
	open(13,file = trim(pdb_id)//"_ATOM.pdb",status = 'new')

	num_atoms = 0
	do
	    read(11,'(a80)',end=10) string
	    if (trim(string(1:4)) == 'ATOM') then
	    	num_atoms = num_atoms + 1
	    	write (13,'(a80)') string
	    end if
	end do
10  close(11)
	close(13)
	
	inquire(file=trim(pdb_id)//"_"//trim(chain_to_read)//"_CA_P_C1_C2.pdb",exist=exist)
	if (exist) then
		open(14,file = trim(pdb_id)//"_"//trim(chain_to_read)//"_CA_P_C1_C2.pdb",status = 'old')
		close(14,status='delete')
	end if
	
	! format:ATOM    110  CD1 PHE A  15      27.005  -7.442  11.099  1.00 19.10           C 
200 format(a4,i7,2x,a3,a1,a3,1x,a1,i4,4x,3f8.3)
	num_chains = len(trim(chain_to_read))
	open(15,file = trim(pdb_id)//"_"//trim(chain_to_read)//"_CA_P_C1_C2.pdb",status = 'new')
	count = 0
	do i = 1, num_chains
		open(16,file = trim(pdb_id)//"_ATOM.pdb",status = 'old')
		do n = 1, num_atoms
			read(16,200) atom_label, atom_num, atom_type, alt_loc, res_type, chain_id, res_num, x, y, z
			if (atom_type(1:2) == 'CA' .or. trim(atom_type) == 'P' .or. &
			trim(atom_type) == "C1'" .or. trim(atom_type) == 'C2') then
				if (alt_loc == ' ' .or. alt_loc == 'A') then
					if (chain_id == chain_to_read(i:i)) then
						count = count + 1
						write(15,200) atom_label, count, atom_type, alt_loc, res_type, chain_id, res_num, x, y, z
					end if
				end if
			end if
		end do
		close(16)
	end do
	close(15)

	open(17,file = trim(pdb_id)//"_ATOM.pdb",status = 'old')
	close(17,status='delete')

end program write_CA_P_C2_C1