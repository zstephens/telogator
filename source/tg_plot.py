import numpy as np
import matplotlib.pyplot as mpl
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

#
#
#
def getColor(i, N, colormap='jet'):
	cm = mpl.get_cmap(colormap) 
	cNorm  = colors.Normalize(vmin=0, vmax=N+1)
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
	colorVal = scalarMap.to_rgba(i)
	return colorVal

#
#
#
def get_read_alignment_polygons(alignments, readlen):
	polygons = []
	p_alpha  = []
	p_color  = []
	p_text   = []
	yp       = [-0.05 * readlen, 0.05 * readlen]
	text_y   = 0.7 * yp[1]
	for i in range(len(alignments)):
		n = alignments[i]
		#
		my_refname = n[2]
		my_refspan = str(n[3])+' - '+str(n[4])
		if my_refname[:3] == 'tel':
			if my_refname[-1] == '?':	# skip these weird regions
				p_color.append('gray')
			else:
				my_refname = 'tel'
				p_color.append('red')
		else:
			p_color.append('blue')
		p_alpha.append(float(n[6]+15.)/(60.+15.))
		#
		xp               = [n[0], n[1]]
		delta_pointy     = 0.8*(n[1]-n[0])
		delta_pointy_rev = 0.2*(n[1]-n[0])
		if n[5] == 'FWD':
			polygons.append(Polygon(np.array([[xp[0],yp[0]], [xp[0],yp[1]], [xp[0]+delta_pointy,yp[1]], [xp[1],0.], [xp[0]+delta_pointy,yp[0]]]), closed=True))
		else:
			polygons.append(Polygon(np.array([[xp[0]+delta_pointy_rev,yp[0]], [xp[0],0.], [xp[0]+delta_pointy_rev,yp[1]], [xp[1],yp[1]], [xp[1],yp[0]]]), closed=True))
		
		p_text.append((n[0], text_y, my_refname+' : '+my_refspan))
		text_y -= 0.40*yp[1]
		#
	axis_val = [0, readlen+1, 2.5*yp[0], 2.5*yp[1]]
	#
	return (polygons, p_color, p_alpha, p_text, axis_val)

#
#
#
def plot_all_read_data(density_data, tl_vals, aln_dat, tel_window, f_title, fig_name, plot_for_paper=False):
	[td_p_e0, td_p_e1, td_q_e0, td_q_e1, p_vs_q_power] = density_data
	[tl_p, tl_q] = tl_vals
	rlen = len(aln_dat[0][7])
	#
	x_off   = tel_window/2
	dens_yt = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
	dens_yl = ['0%', '20%', '40%', '60%', '80%', '100%']
	pred_yt = [-1.0, -0.5, 0.0, 0.5, 1.0]
	#
	if plot_for_paper:
		subplot_nums = [311, 312, 313]
		pred_yl = ['-1', '', '0', '', '1']
		mpl.rcParams.update({'font.size': 20, 'font.weight':'bold'})
		my_figsize = (12,8)
	else:
		subplot_nums = [411, 412, 413]
		pred_yl = ['(q) -1', '', '0', '', '(p) 1']
		mpl.rcParams.update({'font.size': 16, 'font.weight':'bold'})
		my_figsize = (12,10)
	#
	fig = mpl.figure(1, figsize=my_figsize)
	ax1 = mpl.subplot(subplot_nums[0])
	mpl.plot(np.arange(len(td_p_e1))+x_off, td_p_e1, '-m', alpha=0.5)
	mpl.plot(np.arange(len(td_p_e0))+x_off, td_p_e0, '-k', alpha=0.8)
	mpl.rcParams.update({'text.usetex': True})
	mpl.legend(['$p_1(i)$', '$p_0(i)$'], loc=1)
	mpl.rcParams.update({'text.usetex': False})
	mpl.xlim([0,rlen])
	mpl.yticks(dens_yt, dens_yl)
	mpl.ylim([0,1])
	mpl.grid(linestyle='--', alpha=0.5)
	mpl.ylabel('tel density (p)')
	if plot_for_paper == False:
		mpl.title(f_title)
	#
	mpl.subplot(subplot_nums[1], sharex=ax1)
	mpl.plot(np.arange(len(td_q_e1))+x_off, td_q_e1, '-m', alpha=0.5)
	mpl.plot(np.arange(len(td_q_e0))+x_off, td_q_e0, '-k', alpha=0.8)
	mpl.rcParams.update({'text.usetex': True})
	mpl.legend(['$q_1(i)$', '$q_0(i)$'], loc=1)
	mpl.rcParams.update({'text.usetex': False})
	mpl.yticks(dens_yt, dens_yl)
	mpl.ylim([0,1])
	mpl.grid(linestyle='--', alpha=0.5)
	mpl.ylabel('tel density (q)')
	#
	mpl.subplot(subplot_nums[2], sharex=ax1)
	mpl.plot(np.arange(len(p_vs_q_power))+x_off, p_vs_q_power, '-k')
	mpl.rcParams.update({'text.usetex': True})
	mpl.legend(['$S(i)$'], loc=1)
	mpl.rcParams.update({'text.usetex': False})
	mpl.grid(linestyle='--', alpha=0.5)
	mpl.ylabel('tel score')
	mpl.yticks(pred_yt, pred_yl)
	mpl.ylim([-1.0, 1.0])
	polygons = []
	p_color  = []
	p_alpha  = []
	if tl_p > 0:
		xp = [0, tl_p]
		yp = [-1, 1]
		polygons.append(Polygon(np.array([ [xp[0],yp[0]], [xp[0],yp[1]], [xp[1],yp[1]], [xp[1],yp[0]] ]), closed=True))
		p_color.append('red')
		p_alpha.append(0.5)
	if tl_q > 0:
		xp = [rlen - tl_q, rlen]
		yp = [-1, 1]
		polygons.append(Polygon(np.array([ [xp[0],yp[0]], [xp[0],yp[1]], [xp[1],yp[1]], [xp[1],yp[0]] ]), closed=True))
		p_color.append('red')
		p_alpha.append(0.5)
	if len(polygons):
		ax = mpl.gca()
		for i in range(len(polygons)):
			ax.add_collection(PatchCollection([polygons[i]], color=p_color[i], alpha=p_alpha[i]))
	#
	if plot_for_paper == False:
		mpl.subplot(414, sharex=ax1)
		(polygons, p_color, p_alpha, p_text, axis_val) = get_read_alignment_polygons(aln_dat, rlen)
		ax = mpl.gca()
		for i in range(len(polygons)):
			ax.add_collection(PatchCollection([polygons[i]], color=p_color[i], alpha=p_alpha[i]))
		for i in range(len(p_text)):
			mpl.text(p_text[i][0], p_text[i][1], p_text[i][2], ha='left', fontsize=9)
		mpl.axis(axis_val)
		mpl.yticks([],[])
		mpl.grid(linestyle='--', alpha=0.5)
	mpl.xlabel('read position')
	#
	mpl.tight_layout()
	mpl.savefig(fig_name)
	mpl.close(fig)

#
#
#
def tel_len_violin_plot(tel_len_dict, out_fn, plot_means=True, ground_truth_dict={}):
	#
	# plotting constants
	#
	mpl.rcParams.update({'font.size': 18, 'font.weight':'bold'})
	#
	xlab = ['-'] + [str(n) for n in range(1,22+1)] + ['X', 'Y']
	xtck = list(range(1,len(xlab)+1))
	ydel = 5000
	(y_min, y_max) = (20000, 20000)
	ytck  = list(range(-y_min, y_max+ydel, ydel))
	ylab  = []
	for n in ytck:
		if n == 0:
			ylab.append('')
		else:
			ylab.append(str(abs(n)//1000) + 'k')
	#
	ref_2_x = {'chr'+xlab[i]:xtck[i] for i in range(len(xlab))}
	ref_2_x['unanchored'] = xtck[0]
	ref_2_x['unanchore']  = xtck[0]
	#
	v_line_keys     = ['cmeans', 'cmins', 'cmaxes', 'cbars', 'cmedians', 'cquantiles']
	if len(tel_len_dict):
		readcount_denom = max([len(tel_len_dict[k]) for k in tel_len_dict.keys() if k != 'unanchored'])
	else:
		readcount_denom = 1
	width_max       = 1.0
	width_min       = 0.1
	#
	# read in lengths and create data structures needed for violin plot
	#
	(dat_l_p, dat_l_q) = ([], [])
	(dat_p_p, dat_p_q) = ([], [])
	(dat_w_p, dat_w_q) = ([], [])
	for k in tel_len_dict.keys():
		if len(tel_len_dict[k]) == 0:
			continue
		my_width = min( [width_max, max([width_min, width_max*(float(len(tel_len_dict[k]))/readcount_denom)])] )
		if k[-1] == 'p' or k == 'unanchored':
			dat_p_p.append(ref_2_x[k[:-1]])
			dat_l_p.append([])
			dat_w_p.append(my_width)
		elif k[-1] == 'q':
			dat_p_q.append(ref_2_x[k[:-1]])
			dat_l_q.append([])
			dat_w_q.append(my_width)
		for n in tel_len_dict[k]:
			if k[-1] == 'p' or k == 'unanchored':
				dat_l_p[-1].append(n)
			elif k[-1] == 'q':
				dat_l_q[-1].append(-n)
	#
	# violin plot
	#
	fig = mpl.figure(1,figsize=(16,6))
	#
	if len(dat_l_p) and len(dat_p_p):
		violin_parts_p = mpl.violinplot(dat_l_p, dat_p_p, points=200, widths=dat_w_p)
		for pc in violin_parts_p['bodies']:
			pc.set_facecolor('blue')
			pc.set_edgecolor('black')
			pc.set_alpha(0.7)
		for k in v_line_keys:
			if k in violin_parts_p:
				violin_parts_p[k].set_color('black')
				violin_parts_p[k].set_alpha(0.3)
	#
	if len(dat_l_q) and len(dat_p_q):
		violin_parts_q = mpl.violinplot(dat_l_q, dat_p_q, points=200, widths=dat_w_q)
		for pc in violin_parts_q['bodies']:
			pc.set_facecolor('red')
			pc.set_edgecolor('black')
			pc.set_alpha(0.7)
		for k in v_line_keys:
			if k in violin_parts_q:
				violin_parts_q[k].set_color('black')
				violin_parts_q[k].set_alpha(0.3)
	#
	# plot means for each arm
	#
	if plot_means:
		for i in range(len(dat_l_p)):
			yval = np.mean(dat_l_p[i])
			xval = dat_p_p[i]
			mpl.plot([xval - 0.3, xval + 0.3], [yval, yval], '-k', linewidth=2, alpha=0.4)
		for i in range(len(dat_l_q)):
			yval = np.mean(dat_l_q[i])
			xval = dat_p_q[i]
			mpl.plot([xval - 0.3, xval + 0.3], [yval, yval], '-k', linewidth=2, alpha=0.4)
	#
	# plot ground-truth values (for simulated data)
	#
	for k in ground_truth_dict.keys():
		xval = ref_2_x[k[:-1]]
		if k[-1] == 'p':
			yval = ground_truth_dict[k]
		elif k[-1] == 'q':
			yval = -ground_truth_dict[k]
		else:
			print('skipping weird benchmark contig:', k, ground_truth_dict[k])
			continue
		mpl.plot([xval - 0.35, xval + 0.35], [yval, yval], '-k', linewidth=3, alpha=1.0)
	#
	# plot formatting and output
	#
	mpl.plot([0,len(xlab)+1], [0,0], '-k', linewidth=3)
	mpl.xticks(xtck, xlab)
	mpl.xlim([0,len(xlab)+1])
	mpl.yticks(ytck, ylab)
	mpl.ylim([-y_min, y_max])
	mpl.ylabel('<-- q   telomere length   p -->')
	mpl.grid(linestyle='--', alpha=0.5)
	mpl.tight_layout()
	mpl.savefig(out_fn)
	mpl.close(fig)

#
#
#
def anchor_confusion_matrix(conf_dat, out_fn):
	#
	# plotting constants
	mpl.rcParams.update({'font.size': 14, 'font.weight':'bold'})
	#
	xlab = ['-'] + [str(n) for n in range(1,22+1)] + ['X']
	xtck = list(range(1,len(xlab)+1))
	#
	ref_2_ci = {}
	clab     = []
	for ci in xlab:
		if ci == '-':
			ref_2_ci['-'] = 0
			clab.append('-')
		else:
			ref_2_ci['chr' + ci + 'p'] = len(ref_2_ci)
			ref_2_ci['chr' + ci + 'q'] = len(ref_2_ci)
			clab.append('chr' + ci + 'p')
			clab.append('chr' + ci + 'q')
	ref_2_ci['unanchored'] = 0
	#
	Z = [[0 for n in range(len(ref_2_ci))] for m in range(len(ref_2_ci))]
	for k1 in sorted(conf_dat.keys()):
		Z[ref_2_ci[k1[0]]][ref_2_ci[k1[1]]] = conf_dat[k1]
	#
	fig = mpl.figure(3,figsize=(12,10))
	Z = np.array(Z[::-1])
	X, Y = np.meshgrid( range(0,len(Z[0])+1), range(0,len(Z)+1) )
	mpl.pcolormesh(X,Y,Z)
	mpl.axis([0,len(Z[0]),0,len(Z)])
	mpl.yticks(np.arange(0,len(clab))+1.5, clab[::-1])
	mpl.xticks(np.arange(0,len(clab))+0.5, clab, rotation=90)
	mpl.grid(linestyle='--', alpha=0.5)
	mpl.title('subtel confusion matrix')
	mpl.ylabel('ground truth contig')
	mpl.xlabel('where we ended up')
	cb = mpl.colorbar()
	cb.set_label('# reads')
	mpl.tight_layout()
	mpl.savefig(out_fn)
	mpl.close(fig)
