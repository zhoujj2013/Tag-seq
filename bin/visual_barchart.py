import svgwrite
import sys
import os
from svgwrite import cm, mm


boxWidth = 10
box_size = 15
v_spacing = 3

colors = {'G': '#F5F500', 'A': '#FF5454', 'T': '#00D118', 'C': '#26A8FF', 'N': '#B3B3B3'}
cc = {'pp':'#70CDE2', 'pm':'#BE70E2', 'mp': '#E28570','mm':'#94E270'}

def parseDepthFile(infile, strand):
    depth = []
    depth2 = []
    with open(infile, 'r') as f:
        for l in f:
            l_items = l.strip("\n").split('\t')
            d_arr = [l_items[3],l_items[4],l_items[5],l_items[6]]
            depth.append(d_arr)
            pos_sum = int(l_items[3]) + int(l_items[4]) + int(l_items[5]) + int(l_items[6])
            depth2.append(pos_sum)
            #depth2.append(int(l_items[4]))
            #depth2.append(int(l_items[5]))
            #depth2.append(int(l_items[6]))
    if strand == "-":
        depth = depth[::-1]
    
    depth2.sort()
    
    #for i in depth:
    #    print i

    return depth,depth2

def visualizeSites(offtarget_seq, ref_seq, strand, depthFile, prefix):

    # read in depth file
    depth_arr, depth_arr2 = parseDepthFile(depthFile, strand)
    
    # check arrays
    #for i in depth_arr2:
    #    print i
    
    # Initiate canvas
    dwg = svgwrite.Drawing('./' + prefix + '.svg', profile='full')

    if prefix is not None:
        # Define top and left margins
        x_offset = 50
        y_offset = 50
        dwg.add(dwg.text(prefix, insert=(x_offset, 30), style="font-size:20px; font-family:Courier"))
    else:
        # Define top and left margins
        x_offset = 20
        y_offset = 20
        
    x_axis_start_x = x_offset
    x_axis_start_y = y_offset + 180#len(ref_seq) * box_size
    
    x_axis_end_x = x_offset + len(offtarget_seq) * box_size
    x_axis_end_y = y_offset + 180 #len(ref_seq) * box_size

    y_axis_start_x = x_offset
    y_axis_start_y = y_offset

    y_axis_end_x = x_offset
    y_axis_end_y = y_offset + 180 #len(ref_seq) * box_size

    # x-axis
    dwg.add(dwg.line(start=(x_axis_start_x,x_axis_start_y), end=(x_axis_end_x,x_axis_end_y), stroke='black'))

    devide6 = depth_arr2[-1]/6
    devide6_len = len(str(devide6))-1
    base = 10 ** (devide6_len)
    unit = (int(devide6/base) + 1) * base
    
    num_units = (depth_arr2[-1]/unit)+2 # including zero
    
    print devide6
    print base
    print num_units

    for i in range(num_units):
        y_axis_text = i * unit
        y_axis_text_y = y_offset + 182 - i * 30
        y_axis_ticks_y = y_offset + 180 - i * 30
        dwg.add(dwg.text(str(y_axis_text), insert=(x_offset-32, y_axis_text_y), style="font-size:10px; font-family:Courier"))
        dwg.add(dwg.line(start=(x_offset-3,y_axis_ticks_y), end=(x_offset,y_axis_ticks_y), stroke='black'))

    # y-axis
    dwg.add(dwg.line(start=(y_axis_start_x,y_axis_start_y), end=(y_axis_end_x,y_axis_end_y), stroke='black'))

    # Draw bar
    k = 0
    for c in depth_arr:
        pp=c[0]
        pm=c[1]
        mp=c[2]
        mm=c[3]
        x = x_offset + k*box_size
        y = y_offset + 180
        pp_high = (float(pp)/unit)*30
        pm_high = (float(pm)/unit)*30
        mp_high = (float(mp)/unit)*30
        mm_high = (float(mm)/unit)*30

        #print high;

        dwg.add(dwg.rect((x+2, y-pp_high), (box_size-4, pp_high), fill=cc["pp"], stroke="black"))
        dwg.add(dwg.rect((x+2, y-pp_high-pm_high), (box_size-4, pm_high), fill=cc["pm"], stroke="black"))
        dwg.add(dwg.rect((x+2, y-pp_high-pm_high-mp_high), (box_size-4, mp_high), fill=cc["mp"], stroke="black"))
        dwg.add(dwg.rect((x+2, y-pp_high-pm_high-mp_high-mm_high), (box_size-4, mm_high), fill=cc["mm"], stroke="black"))
        k = k + 1

    # Draw legend
    legend_x = x_axis_end_x + 10
    legend_y = y_offset
    dwg.add(dwg.rect((legend_x, legend_y), (box_size, box_size*0.8), fill=cc["pp"], stroke="black"))
    dwg.add(dwg.text('(+) Strand, forward primer', insert=(legend_x + box_size, legend_y+box_size*0.6), style="font-size:12px; font-family:Courier"))
    dwg.add(dwg.rect((legend_x, legend_y + box_size ), (box_size, box_size*0.8), fill=cc["pm"], stroke="black"))
    dwg.add(dwg.text('(+) Strand, reverse primer', insert=(legend_x + box_size, legend_y+box_size+box_size*0.6), style="font-size:12px; font-family:Courier"))
    dwg.add(dwg.rect((legend_x, legend_y + box_size + box_size), (box_size, box_size*0.8), fill=cc["mp"], stroke="black"))
    dwg.add(dwg.text('(-) Strand, forward primer', insert=(legend_x + box_size, legend_y+box_size+box_size+box_size*0.6), style="font-size:12px; font-family:Courier"))
    dwg.add(dwg.rect((legend_x, legend_y + box_size + box_size + box_size), (box_size, box_size*0.8), fill=cc["mm"], stroke="black"))
    dwg.add(dwg.text('(-) Strand, reverse primer', insert=(legend_x + box_size, legend_y+box_size+box_size+box_size+box_size*0.6), style="font-size:12px; font-family:Courier"))

    # Draw x ticks
    y_offset = y_offset + 180 + 11 #len(ref_seq) * box_size + 11
    tick_locations = [1, len(offtarget_seq)]
    tick_locations += range(len(offtarget_seq) + 1)[::10][1:]
    for x in tick_locations:
        dwg.add(dwg.text(str(x), insert=(x_offset + (x - 1) * box_size + 2, y_offset - 2), style="font-size:10px; font-family:Courier"))

    # Draw offtarget row
    for i, c in enumerate(offtarget_seq):
        y = y_offset
        x = x_offset + i * box_size
        dwg.add(dwg.rect((x, y), (box_size, box_size), fill=colors[c]))
        dwg.add(dwg.text(c, insert=(x + 3, y + box_size - 3), fill='black', style="font-size:15px; font-family:Courier"))
    legend_x = x_axis_end_x + 5
    legend_y = y_offset
    dwg.add(dwg.text('Off-target', insert=(legend_x, legend_y+box_size*0.8), style="font-size:15px; font-family:Courier"))
    

    # Draw reference row
    for i, c in enumerate(ref_seq):
        y = y_offset + 30
        x = x_offset + i * box_size
        dwg.add(dwg.rect((x, y), (box_size, box_size), fill=colors[c]))
        dwg.add(dwg.text(c, insert=(x + 3, y + box_size - 3), fill='black', style="font-size:15px; font-family:Courier"))
    legend_x = x_axis_end_x + 5
    legend_y = y_offset + 30
    dwg.add(dwg.text('Reference', insert=(legend_x, legend_y+box_size*0.8), style="font-size:15px; font-family:Courier"))
    
    dwg.save()


def main():
    if len(sys.argv) >= 5:
        # offtarget reference strand depth.bed id
        visualizeSites(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    else:
        print 'Usage: python visualization.py offtarget_seq ref_seq strand count_matrix.txt id'

if __name__ == '__main__':
    main()
