import pyaudi
import random


def get_network(inputs, hidden_layers, units_per_layer):

    weights = []
    bias = []

    prev_layer_outputs = inputs

    #Hidden layers
    for layer in range(hidden_layers):
        this_layer_outputs = []
        layer_weights = []
        layer_bias = []

        for unit in range(units_per_layer):

            unit_output = pyaudi.gdual(0, 'o_{0}_{1}'.format(layer+1, unit), 1)

            b = pyaudi.gdual(1, 'b_{0}_{1}'.format(layer+1, unit), 1)
            unit_output += b

            unit_weights = []
            for prev_idx, prev_output in enumerate(prev_layer_outputs):
                w = pyaudi.gdual(random.normalvariate(0, 1), 'w_{0}_{1}_{2}'.format(layer+1, unit, prev_idx), 1)
                unit_weights.append(w)
                unit_output += w*prev_output

            layer_bias.append(b)
            layer_weights.append(unit_weights)

            unit_output = pyaudi.tanh(unit_output)
            this_layer_outputs.append(unit_output)

        weights.append(layer_weights)
        bias.append(layer_bias)
        prev_layer_outputs = this_layer_outputs

    #output layer
    output = pyaudi.gdual(0, 'N', 1)
    layer = hidden_layers+1
    layer_weights = []
    unit_weights = []
    layer_bias = []

    b = pyaudi.gdual(1, 'b_{0}'.format(layer+1), 1)
    layer_bias.append(b)
    bias.append(layer_bias)
    output += b

    for prev_idx, prev_output in enumerate(prev_layer_outputs):
        w = pyaudi.gdual(random.normalvariate(0, 1), 'w_{0}_{1}'.format(layer, prev_idx), 1)
        unit_weights.append(w)
        output += w*prev_output

    layer_weights.append(unit_weights)
    weights.append(layer_weights)

    return output, weights, bias


def update_weights(NN, w, b, lr):

    dw, db = get_derivatives(NN, w, b)

    for i in range(len(w)):
        for j in range(len(w[i])):
            for k in range(len(w[i][j])):
                w[i][j][k] -= lr*dw[i][j][k]

    for i in range(len(w)):
        for j in range(len(w[i])):
            b[i][j] -= lr*db[i][j]

    return w, b

def update_network(NN, inputs, w, b):

    prev_layer_outputs = inputs
    hidden_layers = len(w)-1
    units_per_layer = len(w[0][0])

    for layer in range(hidden_layers):
        this_layer_outputs = []

        for unit in range(units_per_layer):
            unit_output = pyaudi.gdual(0, 'o_{0}_{1}'.format(layer+1, unit), 1)
            for prev_idx, prev_output in enumerate(prev_layer_outputs):
                unit_output += w[layer][unit][prev_idx] * prev_output
            unit_output = pyaudi.tanh(unit_output)
            this_layer_outputs.append(unit_output)
        prev_layer_outputs = this_layer_outputs

    output = pyaudi.gdual(0, 'N', 1)
    output += b[-1][0]
    for prev_idx, prev_output in enumerate(prev_layer_outputs):
        output += w[layer][0][unit] * prev_output

    return output


def get_derivatives(loss, w, b):

    dw = []
    for wl in w:
        dwl = []
        for wu in wl:
            dwu = []
            for wi in wu:
                dw_idx = [0]*loss.symbol_set_size
                if wi.symbol_set[0] not in loss.symbol_set:
                    dwu.append(0)
                else:
                    idx = loss.symbol_set.index(wi.symbol_set[0])
                    dw_idx[idx] = 1
                    dwu.append(loss.get_derivative(dw_idx))
            dwl.append(dwu)
        dw.append(dwl)

    db = []
    for bl in b:
        dbl = []
        for bu in bl:
            db_idx = [0]*loss.symbol_set_size
            if bu.symbol_set[0] not in loss.symbol_set:
                dbl.append(0)
            else:
                idx = loss.symbol_set.index(bu.symbol_set[0])
                db_idx[idx] = 1
                dbl.append(loss.get_derivative(db_idx))
        db.append(dbl)

    return dw, db


x1 = pyaudi.gdual(0.5, "x1", 1)
x2 = pyaudi.gdual(2, "x2", 1)
x3 = pyaudi.gdual(2, "x2", 1)

y = x1*x2+0.5*x1*x2+pyaudi.sin(x1)*x3

network, weights, biases = get_network([x1, x2], 2, 6)
iters = 10
for i in range(iters):
    loss = ((network - y)**2)
    weights, biases = update_weights(loss, weights, biases, 0.1)
    network = update_network(network, [x1, x2], weights, biases)

print('Loss function value after {0} iterations: {1}\n'.format(iters, ((network - y)**2).constant_cf))
print('y: {0}\nprediction: {1}'.format(y.constant_cf, network.constant_cf))
